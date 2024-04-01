"""


"""

import os
import sys
import fnmatch
import datetime
import logging
import argparse
import yaml
import psutil
from PyQt6 import QtCore
from PyQt6 import QtSql
import CWSaver_splitter


class CWSaver(QtCore.QObject):

    BANNER = """
     _______  _     _  _______  _______  __   __  _______  ______   
    |       || | _ | ||       ||   _   ||  | |  ||       ||    _ |  
    |       || || || ||  _____||  |_|  ||  |_|  ||    ___||   | ||  
    |     __||       || |_____ |       ||       ||   |___ |   |_||_ 
    |    |   |       ||_____  ||       ||       ||    ___||    __  |
    |    |__ |   _   | _____| ||   _   | |     | |   |___ |   |  | |
    |_______||__| |__||_______||__| |__|  |___|  |_______||___|  |_|
                                                                                  
    """

    stop_application = QtCore.pyqtSignal()


    def __init__(self, config_file, parent=None):
        #  initialize the superclass
        super(CWSaver, self).__init__(parent)

        #  store our command line args
        self.config_file = config_file

        self.n_threads = 0
        self.threads = {}
        self.files_to_process = {}
        self.stopping = False
        self.skip_last_file = False

        #  connect the stop signal to our stop method
        self.stop_application.connect(self.stop_app)

        #  create a timer for checking the source directory
        self.file_check_timer = QtCore.QTimer(self)
        self.file_check_timer.setSingleShot(True)
        self.file_check_timer.timeout.connect(self.process_new_files)

        #  start things up after we get the event loop running by using a timer
        QtCore.QTimer.singleShot(0, self.start_app)


    def start_app(self):

        #  add some flair, Joanna
        print(self.BANNER)

        #  read the configuration file
        with open(self.config_file, 'r') as cf_file:
            try:
                self.configuration = yaml.safe_load(cf_file)
            except yaml.YAMLError as exc:
                #  there wa an error reading the config file
                print('Error reading configuration file ' + self.config_file)
                print('  Error string:' + str(exc))
                print('  This application requires a config file.')
                print("Application exiting...")
                QtCore.QCoreApplication.instance().quit()
                return

        #  set up logging
        try:
            #  determine the log level
            if 'log_level' not in self.configuration['application']:
                log_level = logging.INFO
            else:
                log_level = self.configuration['application']['log_level']

            #  get the logging directory, create the dir if needed, and set the log file name
            if 'log_directory' not in self.configuration['application']:
                log_directory = './'
            else:
                log_directory = self.configuration['application']['log_directory']
            log_directory = os.path.normpath(log_directory)
            if not os.path.exists(log_directory):
                os.makedirs(log_directory)
            log_time_string = datetime.datetime.now().strftime("D%Y%m%d-T%H%M%S")
            logfile_name = log_directory + os.sep + 'CWSaver-' + log_time_string + '.log'

            #  create a logger to log to the console
            self.logger = logging.getLogger(__name__)
            self.logger.setLevel(log_level)
            formatter = logging.Formatter('%(asctime)s : %(levelname)s : %(message)s')
            consoleLogger = logging.StreamHandler(sys.stdout)
            fileHandler = logging.FileHandler(logfile_name)
            consoleLogger.setFormatter(formatter)
            fileHandler.setFormatter(formatter)
            self.logger.addHandler(fileHandler)
            self.logger.addHandler(consoleLogger)
            self.logger.info("Staring CWSaver...")
            self.logger.info("Rick Towler - rick.towler@noaa.gov")
            self.logger.info("NOAA/NMFS/AFSC/MACE")

        except Exception as e:
            #  there was an error setting up the log file
            print("Error creating the log file: %s" % (e))
            print("Make sure have write permissions in the logging directory: " + log_directory)
            print("Application exiting...")
            QtCore.QCoreApplication.instance().quit()
            return

        #  open the database file and initialize if needed
        try:
            #  get the path to the sqlite database
            dbFile = (self.configuration['data']['raw_fixed_dest_path'] + os.sep +
                    self.configuration['application']['sqllite_filename'])
            dbFile = os.path.normpath(dbFile)

            self.logger.info("Opening database file: " + dbFile)

            #  check if we're deleting the database file at start
            if self.configuration['application']['purge_database_at_start']:
                self.logger.debug("purge_database_at_start set to True. Deleting " + dbFile)
                try:
                    #  try to delete the file
                    os.remove(dbFile)
                except:
                    #  we're assuming the file doesn't exist
                    pass

            #  open the database file
            self.db = QtSql.QSqlDatabase.addDatabase("QSQLITE", 'CWSaver')
            self.db.setDatabaseName(dbFile)
            if not self.db.open():
                self.logger.critical("Unable to open sqlite database file: " + dbFile)
                self.logger.critical("This file is required. Application exiting...")
                QtCore.QCoreApplication.instance().quit()
                return

            #  check if our tables exist and create if not
            if not 'processed_files' in self.db.tables():
                self.logger.debug("Creating database tables...")
                self.create_database()

        except Exception as e:
            #  there was an error opening or setting up the database
            self.logger.critical("Error opening and initializing the database: %s", e)
            self.logger.exception(e)
            print("Application exiting...")
            QtCore.QCoreApplication.instance().quit()
            return

        #  finish configuring the app and start
        self.logger.info("Configuring processing options...")
        try:

            #  set a bunch of application settings
            self.source_directory = os.path.normpath(self.configuration['data']['raw_file_source_path'])
            self.logger.info("Source file path: " + self.source_directory)

            self.dest_dir = self.configuration['data']['raw_fixed_dest_path']
            self.logger.info("Fixed files will be written to: %s", self.dest_dir)

            self.filename_append = self.configuration['data']['append_to_fixed_filename']
            self.logger.info("append_to_fixed_filename: %s", self.filename_append)

            self.copy_ancillary = self.configuration['data']['copy_ancillary_files']
            self.logger.info("Copy ancillary files: %s", self.copy_ancillary)

            self.max_fixer_threads = self.configuration['application']['max_fixer_threads']
            self.logger.info("Max number of threads to start: %i", self.max_fixer_threads)

            self.overwrite_files = self.configuration['application']['overwrite_files']
            self.logger.info("overwrite_files: %s", self.overwrite_files)

            self.skip_last_file = self.configuration['data']['skip_last_file']
            self.logger.info("skip_last_file: %s", self.skip_last_file)

            self.split_fm = self.configuration['application']['split_fm']
            self.logger.info("split_fm: %s", self.split_fm)

            self.create_instant_echogram = self.configuration['application']['create_instant_echogram']
            self.logger.info("create_instant_echogram: %s", self.create_instant_echogram)
            if self.create_instant_echogram:
                self.instant_echogram_dimensions = [self.configuration['application']['instant_echogram_x'],self.configuration['application']['instant_echogram_y']]

            #  finish setup depending on how we're going to run
            self.monitor = self.configuration['application']['continually_monitor']
            if self.monitor:
                #  the app will continue to run and periodically check and process files in the source dir
                #  when continually monitoring we don't allow reprocessing
                self.reprocess_files = False
                try:
                    # get the directory check interval
                    check_interval = self.configuration['application']['source_check_interval'] * 1000
                except:
                    # the default is 300 seconds
                    check_interval = 300 * 1000

                self.logger.info("continually_monitor: True")
                self.logger.info("reprocess_files: False")
                self.logger.info("The application will continue to run, periodically checking for " +
                        "and processing new files in the source directory.")
                self.logger.info("Press <ctrl>-c to exit.")
                self.logger.info("Setting source directory check interval to " +
                        str(check_interval / 1000) + " seconds.")
                self.file_check_interval = check_interval
                self.file_check_timer.setInterval(500)
                self.file_check_timer.start()
                self.logger.info("Application configured, starting to monitor source directory...")
            else:
                #  only make a single pass thru the source dir and then exit
                self.reprocess_files = self.configuration['application']['reprocess_files']
                self.logger.info("continually_monitor: False")
                self.logger.info("reprocess_files: %s", self.reprocess_files)
                self.logger.info("The application will exit after processing the files in the source directory.")
                self.logger.info("Application configured, Starting a single check of the source directory...")
                self.file_check_timer.setInterval(500)
                self.file_check_timer.start()

        except Exception as e:
            #  there was an error opening or setting up the database
            self.logger.critical("Error configuring processing options: %s", e)
            self.logger.exception(e)
            print("Application exiting...")
            QtCore.QCoreApplication.instance().quit()
            return



    def process_new_files(self):
        '''
        process_new_files will check the source dir for any new files. If found, it will start up to
        max_fixer_threads "fixers" which will read cw/fm file, write the new raw file, and then copy 
        any ancillary files to the destination directory.

        In order to simplify queue and thread management, all new files will be processed before
        this method is called again. In other words, we don't add files to the queue when we are
        currently processing files.

        '''

        #  start a check for new files
        self.process_start_time = datetime.datetime.now()

        #  get the list of raw and ancillary files in the source directory
        self.logger.info("Getting directory listing of source dir: " + self.source_directory)
        new_raw_files = self.get_raw_data_files(self.source_directory)
        scandir_time = datetime.datetime.now() - self.process_start_time
        self.logger.debug("Finished getting directory listing. Elapsed time %s", str(scandir_time))

        #  now filter out the list of files we have already processed
        if not self.reprocess_files:
            #  Get the processed files from the database and remove them from the dict of files to process
            sql = "SELECT file_set_name FROM processed_files"
            query = QtSql.QSqlQuery(sql, self.db)
            while query.next():
                #  get the filename (sans extension)
                this_file = query.value(0)
                #  now try to remove it from the dict of raw files
                new_raw_files.pop(this_file, None)

        #  add any new files to the dict of files to process
        n_new_files = len(new_raw_files)
        self.logger.debug("Found %i new raw files to process", n_new_files)
        self.files_to_process.update(new_raw_files)
        n_files = len(self.files_to_process)

        #  check if we have any files to process and if we can spawn a thread now
        if n_files > 0 and self.n_threads < self.max_fixer_threads:
            self.logger.info("Queue updated. %i file sets in the queue.", n_files)

            #  now spawn fixers, up to max_fixer_threads
            while n_files > 0 and self.n_threads < self.max_fixer_threads:

                #  get a file to process
                basename, file_set = self.files_to_process.popitem()
                n_files -= 1

                #  check if any files in this group are open by another process
                open_file = self.file_is_open(file_set)
                if open_file:
                    #  a file is open and we're not returning open files so we skip
                    self.logger.info("Data file set %s contains an open file: %s. Skipping", basename, open_file)
                    continue

                #  this file is not open so we'll fix it
                self.start_fixer_thread(basename, file_set)

        else:
            #  no files to check, either we wait for more files if we're monitoring or we exit
            if self.monitor:

                #  determine the elapsed processing time
                elapsed_time_ms = (datetime.datetime.now() - self.process_start_time).total_seconds() * 1000
                #  and the number of ms to the next directory check
                next_check_time_ms = int(self.file_check_interval - elapsed_time_ms)
                if next_check_time_ms < 0:
                    next_check_time_ms = 0
                #  nope, no more files, we're monitoring so we start the timer to check the source dir
                self.logger.debug("No new files found. Next check in %i seconds.",
                        (next_check_time_ms / 1000))
                self.file_check_timer.start(next_check_time_ms)

            else:
                #  no new files and we're not monitoring so we'll exit
                self.logger.info("No new files found.")
                self.stop_app()
                return


    def start_fixer_thread(self, basename, file_set):
        '''
        start_fixer_thread will create a CWSaver object and then move it to a
        new thread and start it. The fixer will split out the CW channels from the raw
        file and will optionally copy ancillary files to the destination directory.

        When the fixer finishes, it will signal this application which will
        record the results, retire the fixer/thread and start a new one if there
        are additional files to process.
        '''

        self.logger.info("Starting to fix file set: %s", basename)
        

        #  create a fixer object for the specified file
        fixer = CWSaver_splitter.CWSaver_splitter(basename, file_set, self.dest_dir, self.copy_ancillary,
                self.filename_append, overwrite=self.overwrite_files,split_fm=self.split_fm, 
                create_instant_echogram=self.create_instant_echogram, instant_echogram_dimensions=self.instant_echogram_dimensions)

        #  create a thread to run the fixer in
        thread = QtCore.QThread()

        #  move the fixer to it
        fixer.moveToThread(thread)

        #  connect thread specific signals and slots - this facilitates starting,
        #  stopping, and deletion of the threads.
        thread.started.connect(fixer.start_processing)
        fixer.finished.connect(self.fixer_finished)
        fixer.file_error.connect(self.fixer_error)
        fixer.file_read_progress.connect(self.fixer_progress)
        fixer.file_write_progress.connect(self.fixer_progress)
        thread.finished.connect(thread.deleteLater)

        #  increment our thread counter
        self.n_threads += 1

        #  store references to our new objects
        self.threads[fixer] = thread

        #  and finally, start the thread - this will also start polling
        thread.start()


    @QtCore.pyqtSlot(object,str, str,datetime.datetime, str, int, int, int, bool, int, object, int)
    def fixer_finished(self, fixer_obj, basename, in_file, processing_time, out_file,
            bytes_read, bytes_written, n_channels, split_fm, file_error, elapsed_seconds, n_ancillary_copied):

        self.logger.info("File set split: %s  :: %i CW Channels ::  elapsed time %s", basename,
                n_channels, str(elapsed_seconds))

        #  write the results to the database
        time_str = self.datetime_to_db_str(processing_time)
        sql = ("INSERT INTO processed_files VALUES('" + basename + "','" + time_str + "','" + in_file + "','" + out_file + "'," +
                str(bytes_read) + "," + str(bytes_written) + "," + str(n_channels) + "," + str(split_fm) + "," +str(file_error) + ",'" +
                str(elapsed_seconds) + "'," + str(n_ancillary_copied) + ")")
        query = QtSql.QSqlQuery(sql, self.db)
        query.exec()

        #  clean up this thread and decide the next move
        self.cycle_fixer(fixer_obj)


    @QtCore.pyqtSlot(object, str, str, datetime.datetime, object)
    def fixer_error(self, fixer_obj, fileset_name, in_file, processing_time, error_obj):

        #  log the error
        self.logger.error("Error processing file: %s", in_file)
        self.logger.error("Error string: %s", str(error_obj))

        #  write the results to the database
        time_str = self.datetime_to_db_str(processing_time)
        sql = ("INSERT INTO errors VALUES('" + fileset_name + "','" + time_str + "','" + in_file + "','" + str(error_obj) + "')")
        query = QtSql.QSqlQuery(sql, self.db)
        query.exec()


    def cycle_fixer(self, fixer_obj):
        '''
        cycle_fixer will clean up after a fixer finishes or errors out and if
        there are files to process (and we have not been told to stop) it will
        start a new fixer to fix the next available file. If there are no files
        to fix it will start the timer to
        '''

        #  get a reference to this fixer's thread, stop it, and remove the fixer from
        #  the thread mapping dict,
        thread = self.threads[fixer_obj]
        thread.quit()
        del self.threads[fixer_obj]
        self.n_threads -= 1
        self.logger.debug("Current thread count: %i", self.n_threads)

        #  determine our next move
        if not self.stopping:

            #  we've not been told to stop so check if we have any more files to process
            n_files = len(self.files_to_process)
            if n_files > 0:
                #  yes, continue working thru them until we find one to process
                for i in range(n_files):

                    #  get a file to process
                    basename, file_set = self.files_to_process.popitem()

                    #  check if any files in this group are open by another process
                    open_file = self.file_is_open(file_set)
                    if open_file:
                        #  a file is open and we're not returning open files so we skip
                        self.logger.info("Data file set %s contains an open file: %s. Skipping", basename, open_file)
                        continue

                    #  this file is not open so we'll fix it
                    self.start_fixer_thread(basename, file_set)

                    #  we're done here
                    return

            else:
                if self.monitor:
                    if self.n_threads == 0:
                        #  determine the elapsed processing time
                        elapsed_time_ms = (datetime.datetime.now() - self.process_start_time).total_seconds() * 1000
                        #  and the number of ms to the next directory check
                        next_check_time_ms = int(self.file_check_interval - elapsed_time_ms)
                        if next_check_time_ms < 0:
                            next_check_time_ms = 0
                        #  nope, no more files, we're monitoring so we start the timer to check the source dir
                        self.logger.debug("All files in the queue have been processed. Next check in %i seconds.",
                                (next_check_time_ms / 1000))
                        self.file_check_timer.start(next_check_time_ms)
                else:
                    #  we're not monitoring and we're out of files - check if we can exit
                    if self.n_threads == 0:
                        #  yes, we can stop
                        self.logger.info("All files have been processed.")
                        self.stop_app(finished=True)
                    else:
                        #  nope, we're waiting for some other threads to finish
                        self.logger.debug("No more files to process but waiting for %i thread(s) to finish", self.n_threads)

        else:
            #  we've been told to stop but we need to wait for all threads to finish.
            #  if n_threads == 0, we will call stop_app to exit
            if self.n_threads == 0:
                self.stop_app(finished=True)
            else:
                self.logger.debug("Waiting for %i thread(s) to finish before exiting", self.n_threads)



    @QtCore.pyqtSlot(str, int)
    def fixer_progress(self, in_file, cumulative_pct):
        pass



    def get_raw_data_files(self, path, return_bot=False, require_idx=True):
        '''
        get_raw_data_files returns a dict of dicts, keyed by raw file base name, of raw
        files and their associated data files. Each value contains a dict containing the
        full path to the raw file, the idx file, all .xyz files, and optionally
        the bot file. For convenience, it also contains a list of all of the ancillary
        files, making it easier to operate on the files as a group.
        '''

        #  initialize the return dictionary
        raw_data_fileset = {}

        #  use scandir to get a directory listing
        files = os.scandir(path)

        #  create a list of files from the scandir output
        files = [os.path.normpath(file.path) for file in files]

        #  get the raw files
        raw_files = fnmatch.filter(files, '*.raw')

        #  and sort
        raw_files.sort()

        #  determine how many raw files we have
        n_raw_files = len(raw_files)

        #  if we're skipping the last file, reduce the count by one
        if self.skip_last_file:
            n_raw_files -= 1

        #  iterate thru the raw files extracting info about the ancillary files and
        #  build the return dict
        for n in range(n_raw_files):
            #  get the full file name
            fullname = raw_files[n]

            #  split the file base name and extension
            filepath, filename = os.path.split(fullname)
            basename, ext = os.path.splitext(filename)

            #  now get the ancillary files for this raw file
            xyz_files = fnmatch.filter(files, '*' + basename + '*.XYZ')
            idx_file = fnmatch.filter(files, '*' + basename + '.idx')
            if len(idx_file) > 0:
                idx_file = idx_file[0]
            else:
                idx_file = None
            if return_bot:
                bot_file = fnmatch.filter(files, '*' + basename + '.bot')
                if len(bot_file) > 0:
                    bot_file = bot_file[0]
                else:
                    bot_file = None

            #  check for .idx file if required
            if require_idx and len(idx_file) == 0:
                #  no idx file, skip this file
                self.logger.debug("Data file %s is missing the .idx file. Skipping", basename)
                continue

            #  set the entry for this file in our return dict
            raw_data_fileset[basename] = {}
            raw_data_fileset[basename]['.raw'] = fullname
            raw_data_fileset[basename]['.idx'] = idx_file
            raw_data_fileset[basename]['.xyz'] = xyz_files
            raw_data_fileset[basename]['ancillary_list'] = [idx_file]
            raw_data_fileset[basename]['ancillary_list'].extend(xyz_files)
            if return_bot:
                raw_data_fileset[basename]['ancillary_list'].append(bot_file)

        return raw_data_fileset


    def file_is_open(self, file_set):
        '''
        Returns true if another process any files from the file set open

        Credit to Stack Overflow user Tavy
        https://stackoverflow.com/questions/11114492/check-if-a-file-is-not-open-nor-being-used-by-another-process
        '''

        for proc in psutil.process_iter():
            try:
                for item in proc.open_files():
                    if file_set['.raw'] == item.path:
                        return file_set['.raw']
                    for f in file_set['ancillary_list']:
                        if f == item.path:
                            return f
            except Exception:
                pass

        return ''


    def create_database(self):
        '''
        create_database creates the tables in the output tracking database
        '''

        # list of SQL statements that define the base camtrawlMetadata database schema
        sql = ["CREATE TABLE processed_files (file_set_name TEXT NOT NULL, time TEXT NOT NULL, source_file TEXT, " +
                "fixed_file TEXT, bytes_read INTEGER, bytes_written INTEGER, n_channels INTEGER , split_fm BOOLEAN , error INTEGER, " +
                "processing_time TEXT, n_ancillary_files_copied INTEGER, PRIMARY KEY(file_set_name,time))",
               "CREATE TABLE errors (file_set_name TEXT NOT NULL, time TEXT NOT NULL, source_file TEXT NOT NULL, " +
                    "error TEXT NOT NULL, PRIMARY KEY(file_set_name,time))"]

        #  execute the sql statements
        for s in sql:
            query = QtSql.QSqlQuery(s, self.db)
            query.exec()


    def datetime_to_db_str(self, dt_obj):

        dt_string = dt_obj.strftime("%Y-%m-%d %H:%M:%S")
        dt_string = dt_string + '.%03d' % (round(dt_obj.microsecond / 1000))

        return dt_string


    @QtCore.pyqtSlot()
    def stop_app(self, finished=False):

        if not finished:
            self.logger.info("Starting to shut down CWSaver...")
            self.file_check_timer.stop()

        if self.n_threads > 0:
            self.logger.info("Waiting for fixer threads to finish. This may take a bit...")
            self.stopping = True
        else:
            #  we can stop immediately since there are no fixer threads running
            if self.db.isOpen():
                self.logger.debug("Closing database...")
                self.db.close()

            self.logger.info("Application exiting...")
            QtCore.QCoreApplication.instance().quit()
            return


    def external_stop(self):
        '''
        external_stop is called when one of the main thread exit handlers are called.
        It emits a stop signal that is then received by the console_app which then
        shuts everything down in the application's thread.
        '''
        self.stop_application.emit()


def exit_handler(a,b=None):
    '''
    exit_handler is called when CTRL-c is pressed on Windows
    '''
    global ctrlc_pressed

    if not ctrlc_pressed:
        #  make sure we only act on the first ctrl-c press
        ctrlc_pressed = True
        print("CTRL-C detected. Shutting down...")
        console_app.external_stop()

    return True


def signal_handler(*args):
    '''
    signal_handler is called when ctrl-c is pressed when the python console
    has focus. On Linux this is also called when the terminal window is closed
    or when the Python process gets the SIGTERM signal.
    '''
    global ctrlc_pressed

    if not ctrlc_pressed:
        #  make sure we only act on the first ctrl-c press
        ctrlc_pressed = True
        print("CTRL-C or SIGTERM/SIGHUP detected. Shutting down...")
        console_app.external_stop()

    return True


if __name__ == '__main__':

    #  create a state variable to track if the user typed ctrl-c to exit
    ctrlc_pressed = False

    config_file = './CWSaver.yml'

    #  Set up the handlers to trap ctrl-c
    if sys.platform == "win32":
        #  On Windows, we use win32api.SetConsoleCtrlHandler to catch ctrl-c
        import win32api
        win32api.SetConsoleCtrlHandler(exit_handler, True)
    else:
        #  On linux we can use signal to get not only ctrl-c, but
        #  termination and hangup signals also.
        import signal
        signal.signal(signal.SIGINT, signal_handler)
        signal.signal(signal.SIGTERM, signal_handler)
        signal.signal(signal.SIGHUP, signal_handler)

    #  parse the command line arguments
    parser = argparse.ArgumentParser(description='CWSaver creates new raw files containing only CW channels.')
    parser.add_argument("-f", "--config_file", help="Specify the path to the yml configuration file.")

    args = parser.parse_args()

    if (args.config_file):
        config_file = os.path.normpath(str(args.config_file))

    #  create an instance of QCoreApplication and an instance of the our application
    app = QtCore.QCoreApplication(sys.argv)
    console_app = CWSaver(config_file, parent=app)

    #  and start the event loop
    sys.exit(app.exec())
