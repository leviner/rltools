import numpy as np


class Calculation():
    #
    # Chapter IIB: Signal generation
    #

    @staticmethod
    def generateIdealWindowedTransmitSignal(f0, f1, tau, fs, slope):
        """
        Generate the ideal windowed transmit signal.
        
        Parameters
        ----------
        f0 : float
            The start frequency [Hz]
        f1 : float
            The end frequency [Hz]
        tau : float
            The transmit signal duration [s]
        fs : float
            The signal sample rate [Hz]
        slope : float
            The proportion of the pulse that is windowed [1]. This includes the start and end parts of the windowing.
                 
        Returns
        -------
        y : np.array
            The ideal windowed transmit signal amplitude [1]
        t : np.arrayfloat
            The time of each value in y [s]
    
        """
        
        nsamples = int(np.floor(tau * fs))
        t = np.linspace(0, nsamples - 1, num=nsamples) * 1 / fs
        y = Calculation.chirp(t, f0, tau, f1)
        L = int(np.round(tau * fs * slope * 2.0))  # Length of hanning window
        w = Calculation.hann(L)
        N = len(y)
        w1 = w[0 : int(len(w) / 2)]
        w2 = w[int(len(w) / 2) :]
        i0 = 0
        i1 = len(w1)
        i2 = N - len(w2)
        i3 = N
        y[i0:i1] = y[i0:i1] * w1
        y[i2:i3] = y[i2:i3] * w2

        return y, t

    @staticmethod
    def chirp(t, f0, t1, f1):
        """
        Generate a chirp pulse.
        
        Parameters
        ----------
        t : float
            Times at which to calculate the chirp signal [s]
        f0 : float
            Start frequency of the pulse [Hz]
        t1 : float
            Duration of the pulse [s]
        f1 : float
            End frequency of the pulse [Hz]
            
        Returns
        -------
        np.array
            The chirp pulse amplitude at times t [1]
        """
        
        a = np.pi * (f1 - f0) / t1
        b = 2 * np.pi * f0

        return np.cos(a * t * t + b * t)

    @staticmethod
    def hann(L):
        """
        Generate Hann window weights.
        
        Parameters
        ----------
        L : int
            The number of samples to use [1]
        
        Returns
        -------
        np.array
            The Hann window weights [1]
        """
        
        n = np.arange(0, L, 1)

        return 0.5 * (1.0 - np.cos(2.0 * np.pi * n / (L - 1)))

    #
    # Chapter IIC: Signal reception
    #

    @staticmethod
    def calcDecmiatedSamplingRate(filter_v, f_s):
        """
        Calculate the decimated sampling frequency after each filter stage.
        
        Parameters
        ----------
        filter_v : sequence of np.arrays
            The filters.
        f_s : float
            The undecimated sample rate [Hz]
            
        Returns
        -------
        np.array
            Decimated sampled rate after each given filter stage [Hz]
        """
        
        f_s_dec = [f_s]
        if filter_v is not None:
            for v, _filter_v in enumerate(filter_v):
                f_s_dec.append(f_s_dec[v] / _filter_v["D"])

        return f_s_dec

    #
    # Chapter IID: Pulse compression
    #

    @staticmethod
    def calcNormalizedTransmitSignal(y_tx_n):
        """
        Normalise the transmit signal by the maximum of the transmit signal.
        
        Parameters
        ----------
        y_tx_n : np.array
            The transmit signal [V]
        
        Returns
        -------
        np.array
            The normalised transmit signal [1]
            
        """
        
        return y_tx_n / np.max(y_tx_n)

    @staticmethod
    def calcFilteredAndDecimatedSignal(y_tilde_tx_n, filter_v):
        """
        Filter and decimate a signal using given filters.
        
        Parameters
        ----------
        y_tilde_tx_n : np.array
            Normalised transmit signal [1]
        filter_v : sequence of np.arrays
            The filter coefficients for each filter to apply [1]
            
        Returns
        -------
        np.array
            The filtered and decimated transmit signal [1]
    """
        
        # Initialize with normalized transmit pulse
        y_tilde_tx_nv = [y_tilde_tx_n]
        v = 0

        if filter_v is not None:
            for filter_vi in filter_v:
                tmp = np.convolve(y_tilde_tx_nv[v], filter_vi["h_fl_i"], mode="full")[
                    0 :: filter_vi["D"]
                ]
                y_tilde_tx_nv.append(tmp)
                v += 1

        return y_tilde_tx_nv

    @staticmethod
    def calcAutoCorrelation(y_mf_n, f_s_dec):
        """
        Calculate the autocorrelation of the matched filter signal.
        
        Parameters
        ----------
        y_mf_n : np.array
            Input signal [1]
        f_s_dec : float
            Decimated sampling frequency [Hz]
            
        Returns
        -------
        y_mf_auto_n : np.array
            Autocorrelation of the input filter [1]
        t_eff : float
            Effective transmit pulse duration [s]
        """
        
        y_mf_n_conj_rev = np.conj(y_mf_n)[::-1]
        y_mf_twoNormSquared = np.linalg.norm(y_mf_n, 2) ** 2
        y_mf_n_conj_rev = y_mf_n_conj_rev
        y_mf_twoNormSquared = y_mf_twoNormSquared

        # Calculate auto correlation function for matched filter
        y_mf_auto_n = np.convolve(y_mf_n, y_mf_n_conj_rev) / y_mf_twoNormSquared

        # Estimate effective pulse duration
        p_tx_auto = np.abs(y_mf_auto_n) ** 2
        tau_eff = np.sum(p_tx_auto) / ((np.max(p_tx_auto)) * f_s_dec)

        return y_mf_auto_n, tau_eff

    @staticmethod
    def calcPulseCompressedSignals(quadrant_signals, y_mf_n):
        """
        Pulse compress input signals.
        
        Parameters
        ----------
        quadrant_signals : sequence of np.arrays
            The signals to be pulse compressed [V]
        y_mf_n : np.array
            The matched filter to apply to the input signals [1]
            
        Returns
        -------
        np.array
           The pulse compressed signals [V]
        """
        
        # Do pulse compression on all quadrants
        y_mf_n_conj_rev = np.conj(y_mf_n)[::-1]
        y_mf_twoNormSquared = np.linalg.norm(y_mf_n, 2) ** 2

        pulseCompressedQuadrants = []
        start_idx = len(y_mf_n_conj_rev) - 1
        for u in quadrant_signals:
            # Pulse compression per quadrant
            y_pc_nu = np.convolve(y_mf_n_conj_rev, u, mode="full") / y_mf_twoNormSquared
            # Correct sample indexes for mached filter
            y_pc_nu = y_pc_nu[start_idx::]
            pulseCompressedQuadrants.append(y_pc_nu)

        return np.array(pulseCompressedQuadrants)

    @staticmethod
    def calcAverageSignal(y_pc_nu):
        """
        Calculate the average pulse compressed signal over all 
        receivers/transducer sectors.
        
        Parameters
        ----------
        y_pc_nu : np.array
            The pulse compressed data from each receiver/transducer sector [V]
        Returns
        -------
        np.array
            The average pulse compressed signal [V]
        """
        
        return np.sum(y_pc_nu, axis=0) / y_pc_nu.shape[0]

    @staticmethod
    def calcTransducerHalves(y_pc_nu):
        """
        Calculate the half transducer pulse compressed signals from the four 
        sectors in a 4-sector transducer.
        
        Parameters
        ----------
        y_pc_nu : np.array
            The pulse compressed data from each receiver/transducer sector [V]
        
        Returns
        -------
        y_pc_fore_n : np.array
            Signal from the forward half of the transducer [V]
        y_pc_aft_n : np.array
            Signal from the aft half of the transducer [V]
        y_pc_star_n : np.array
            Signal from the starboard half of the transducer [V]
        y_pc_port_n : np.array
            Signal from the port half of the transducer [V]
        """
        y_pc_fore_n = 0.5 * (y_pc_nu[2, :] + y_pc_nu[3, :])
        y_pc_aft_n = 0.5 * (y_pc_nu[0, :] + y_pc_nu[1, :])
        y_pc_star_n = 0.5 * (y_pc_nu[0, :] + y_pc_nu[3, :])
        y_pc_port_n = 0.5 * (y_pc_nu[1, :] + y_pc_nu[2, :])

        return y_pc_fore_n, y_pc_aft_n, y_pc_star_n, y_pc_port_n

    #
    # Chapter IIE: Power angles and samples
    #

    @staticmethod
    def calcPower(y_pc, z_td_e, z_rx_e, nu):
        """
        Calculate the received power into a matched load.
        
        Output received power values of 0.0 are set to 1e-20.
        
        Parameters
        ----------
        y_pc : np.array
            Pulse compressed signal [V]
        z_td_e : float
            Transducer electrical impedance [Ω]
        z_rx_e : float
            Receiver electrical impedance [Ω]
        nu : int
            Number of receiver channels [1]
        
        Returns
        -------
        np.array
            Received electrical power [W]
            
        """
        K1 = nu / ((2 * np.sqrt(2)) ** 2)
        K2 = (np.abs(z_rx_e + z_td_e) / z_rx_e) ** 2
        K3 = 1.0 / np.abs(z_td_e)
        C1Prx = K1 * K2 * K3
        Prx = C1Prx * np.abs(y_pc) ** 2
        Prx[Prx == 0] = 1e-20

        return Prx

    @staticmethod
    def calcAngles(y_pc_halves, gamma_theta, gamma_phi):
        """
        Calculate splitbeam angles from 4 quadrant receiver/transducer data.
        
        Parameters
        ----------
        y_pc_halves : np.array
            Pulse compressed signal from transducer halves [V]
        gamma_theta : float
             Major axis conversion factor from electrical splitbeam angles to physical angles []
        gamma_phi : float
            Minor axis conversion factor from electrical splitbeam angles to physical angles []
        
        Returns
        -------
        theta_n : np.array
            Major axis physical splitbeam angles [°]
        phi_n : np.array
            Minor axis physical splitbeam angles [°]
        
        """

        y_pc_fore_n, y_pc_aft_n, y_pc_star_n, y_pc_port_n = y_pc_halves

        y_theta_n = y_pc_fore_n * np.conj(y_pc_aft_n)
        y_phi_n = y_pc_star_n * np.conj(y_pc_port_n)

        theta_n = (
            np.arcsin(np.arctan2(np.imag(y_theta_n), np.real(y_theta_n)) / gamma_theta)
            * 180
            / np.pi
        )
        phi_n = (
            np.arcsin(np.arctan2(np.imag(y_phi_n), np.real(y_phi_n)) / gamma_phi)
            * 180
            / np.pi
        )

        return theta_n, phi_n

    #
    # Chapter III: TARGET STRENGTH
    #

    @staticmethod
    def calcSp(p_rx_e_n, r_n, alpha_f_c, p_tx_e, lambda_f_c, g_0_f_c, r0=None, r1=None):
        """
        Calculate the point scattering strength.
        
        Parameters
        ----------
        p_rx_e_n : np.array
            Received electric power into a matched load [W]
        r_n : np.array
            Range of samples in `p_rx_e_n` [m]
        alpha_f_c : float
            Acoustic absorption at centre frequency of transmitted chirp signal [Hz]
        lambda_f_c : float
            Acoustic wavelength at centre frequency of transmitted chirp signal [m]
        g_0_f_c : 
            Transducer gain centre frequency of transmitted chirp signal [1]
        r0 : float, optional
            Start range (inclusive) of data to be processed [m]
        r1 : float, optional
            End range of (inclusive) data to be processed [m]
            
        Returns
        -------
        np.array
            Point scattering strength [dB re 1m^2]
        """
        
        # Pick a range of the data
        if r0 is not None and r1 is not None:
            Idx = np.where((r_n >= r0) & (r_n <= r1))
            r_n = r_n[Idx]
            p_rx_e_n = p_rx_e_n[Idx]

        # Calculate Sp for the choosen range
        S_p_n = (
            10.0 * np.log10(p_rx_e_n)
            + 40.0 * np.log10(r_n)
            + 2.0 * alpha_f_c * r_n
            - 10
            * np.log10((p_tx_e * lambda_f_c**2 * g_0_f_c**2) / (16.0 * np.pi**2))
        )

        return S_p_n

    @staticmethod
    def singleTarget(
        y_pc_n, p_rx_e_n, theta_n, phi_n, r_n, r0, r1, before=0.5, after=0.5
    ):
        """
        Detect single targets in acoustic data.
        
        This is a pseudo single target detector (SED) using the max peak
        power within the range interval [r0 r1] as the single target
        detection criteria.

        Parameters
        ----------
        y_pc_n : np.array
            Average pulse compressed signal over all sectors [V]
        p_rx_e_n : np.array
            Received electric power into a matched load [W]
        theta_n : np.array
            Minor axis splitbeam angles of data in `p_rx_e_n` [°]
        phi_n : np.array
            Major axis splitbeam angles of data in `p_rx_e_n` [°]
        r_n :  np.array
            Range of samples in `p_rx_e_n` [m]
        r0 : float
            Start (inclusive) of range interval for SED [m]
        r1 : float
            End (inclusive) of range interval for SED [m]
        before : float, optional
            Distance before the detected target to extract [m]
        after : float, optional
            Distance after the detected target to extract [m]
            
        Returns
        -------
        r_t : float
            Range to peak of detected target [m]
        theta_t : float
            Minor axis splitbeam of detected target [°]
        phi_t : float
            Major axis splitbeam of detected target [°]
        y_pc_t : np.array
            Pulse compressed signal from `before` to `after' of detected target [V]
        p_rx_t : np.array
            Received electric power from `before` to `after' of detected target [W]
        dum_theta : np.array
            Minor axis splitbeam angles from `before` to `after' of detected target [°]
        dum_phi : np.array
            Major axis splitbeam angles from `before` to `after' of detected target [°]
        dum_r : np.array
            Range from `before` to `after' of detected target [m]
        """

        # Get the samples within the sub range defined by the range interval
        if r0 is not None and r1 is not None:
            Idx = np.where((r_n >= r0) & (r_n <= r1))
            r_n_sub = r_n[Idx]
            y_pc_n_sub = y_pc_n[Idx]
            p_rx_e_n_sub = p_rx_e_n[Idx]
            theta_n_sub = theta_n[Idx]
            phi_n_sub = phi_n[Idx]

        # Use peak power within given range limits as index for single target
        # to get the range, theta and phi for the single target
        idx_peak_p_rx = np.argmax(p_rx_e_n_sub)
        r_t = r_n_sub[idx_peak_p_rx]
        theta_t = theta_n_sub[idx_peak_p_rx]
        phi_t = phi_n_sub[idx_peak_p_rx]

        # Extract pulse compressed samples "before" and "after" the peak power
        r_t_begin = r_t - before
        r_t_end = r_t + after
        Idx2 = np.where((r_n_sub >= r_t_begin) & (r_n_sub <= r_t_end))
        y_pc_t = y_pc_n_sub[Idx2]
        p_rx_t = p_rx_e_n_sub[Idx2]
        dum_theta = theta_n_sub[Idx2]
        dum_phi = phi_n_sub[Idx2]
        dum_r = r_n_sub[Idx2]

        return r_t, theta_t, phi_t, y_pc_t, p_rx_t, dum_theta, dum_phi, dum_r

    @staticmethod
    def alignAuto(y_mf_auto_n, y_pc_t_n):
        """
        Calculate the reduced autocorrelation of a matched filter.
        
        Parameters
        ----------
        y_mf_auto_n : np.array
            Autocorrelation of the matched filter [1]
        y_pc_t_n : np.array
            Pulse compressed signal from a single target [V]
            
        Returns
        -------
        np.array
            Reduced autocorrelation of the matched filter [1]
        
        """
        # The equivalent samples around the peak in the decimated tx signal
        idx_peak_auto = np.argmax(np.abs(y_mf_auto_n))
        idx_peak_y_pc_t_n = np.argmax(np.abs(y_pc_t_n))

        left_samples = idx_peak_y_pc_t_n
        right_samples = len(y_pc_t_n) - idx_peak_y_pc_t_n

        idx_start_auto = max(0, idx_peak_auto - left_samples)
        idx_stop_auto = min(len(y_mf_auto_n), idx_peak_auto + right_samples)

        y_mf_auto_red_n = y_mf_auto_n[idx_start_auto:idx_stop_auto]

        return y_mf_auto_red_n

    @staticmethod
    def calcDFTforTS(y_pc_t_n, y_mf_auto_red_n, n_f_points, f_m, f_s_dec):
        """
        Calculate the discrete Fourier transforms (DFT).
        
        Provides the DFT of the target signal, DFT of the reduced auto 
        correlation signal, and the normalized DFT of the target signal.
        
        Parameters
        ----------
        y_pc_t_n : 
        y_mf_auto_red_n :
        n_f_points :
            Maximum number of points to use in the DFT (actual number will be
            the largest power of two that is less than n_f_points) [1]
        f_m : np.array
            Frequency vcetor [Hz]
        f_s_dec : float
            Decimated sampling rate [Hz]
        
        Returns
        -------
        Y_pc_t_m :
            DFT of the pulse compressed signal from a single target [V]
        Y_mf_auto_red_m :
            DFT of the reduced autocorrelation function for the matched filter [1]
        Y_tilde_pc_t_m :
            DFT of `Y_pc_t_m` normalised by `Y_mf_auto_red_m` [1]
        
        """
        # The number of DFT points inpower of 2
        N_DFT = int(2 ** np.ceil(np.log2(n_f_points)))

        idxtmp = np.floor(f_m / f_s_dec * N_DFT).astype("int")
        idx = np.mod(idxtmp, N_DFT)

        # DFT for the target signal
        _Y_pc_t_m = np.fft.fft(y_pc_t_n, n=N_DFT)
        Y_pc_t_m = _Y_pc_t_m[idx]

        # DFT for the transmit signal
        _Y_mf_auto_red_m = np.fft.fft(y_mf_auto_red_n, n=N_DFT)
        Y_mf_auto_red_m = _Y_mf_auto_red_m[idx]

        # The Normalized DFT
        Y_tilde_pc_t_m = Y_pc_t_m / Y_mf_auto_red_m

        return Y_pc_t_m, Y_mf_auto_red_m, Y_tilde_pc_t_m

    @staticmethod
    def calcPowerFreqTS(N_u, Y_tilde_pc_t_m, z_td_e, z_rx_e):
        """
        Calculate the received power spectrum for a single target
        
        Parameters
        ----------
        N_u : int
            Number of transducer sectors/receiver channels
        Y_tilde_pc_t_m :
            DFT of the pulse compressed signal from a single target normalised 
            by the DFT of the reduced autocorrelation function of the matched
            filter [1]
        z_td_e :
            Transducer sector electric impedance [Ω]
        z_rx_e :
            Receiver electric impedance [Ω]
        
        Returns
        -------
        np.array
            DFT of received electric power in a matched load for the signal from 
            a single target
            
        """
        imp = (np.abs(z_rx_e + z_td_e) / np.abs(z_rx_e)) ** 2 / np.abs(z_td_e)
        P_rx_e_t_m = N_u * (np.abs(Y_tilde_pc_t_m) / (2 * np.sqrt(2))) ** 2 * imp
        return P_rx_e_t_m

    @staticmethod
    def calcg0_calibrated(f, freq, gain):
        """
        Interpate frequency dependent gain data onto specified frequencies.
        
        Uses one-dimensional piecewise linear interpolation.
        
        Parameters
        ----------
        f : np.array
            Frequencies to interpolate on to [Hz]
        freq : np.array
            Frequencies at which gain values are provided [Hz]
        gain : np.array
            Gain values to use in the interpolation [dB]
        
        Returns
        -------
        np.array
            The gain values at frequencies given in `f` [1]
        
        """
        dB_G0 = np.interp(f, freq, gain)

        return np.power(10, dB_G0 / 10)

    @staticmethod
    def calcg0_notcalibrated(f, f_n, G_f_n):
        """
        Estimate frequency dependent gain based on theoretical variation with frequency.
        
        Parameters
        ----------
        f : np.array
            Frequencies at which to calculate the gain [Hz]
        f_n : float
            Frequency that the given gain is for [Hz]
        G_f_n : float
            Gain at `f_n` [dB]
        
        Returns
        -------
        np.array
            The gain values at frequencies given in `f` [1]
        """
        dB_G0 = G_f_n + 20 * np.log10(f / f_n)

        return np.power(10, dB_G0 / 10)

    @staticmethod
    def calcTSf(P_rx_e_t_m, r_t, alpha_m, p_tx_e, lambda_m, g_theta_t_phi_t_f_t):
        """
        Calculate TS as a function of frequency for a single target.
        
        Parameters
        ----------
        P_rx_e_t_m : np.array
            DFT of the received electric power in a matched load for the signal
            from a single target at a range of frequencies [W]
        r_t : float
            Range to the target [m]
        alpha_m : np.array
            Acoustic absorption at a range of frequencies [dB/m]
        p_tx_e : float
            Transmitted electric power [W]
        lambda_m : 
            Acoustic wavelength at a range of frequencies [m]
        g_theta_t_phi_t_f_t : np.array
            Transducer gain as a function of echo arrival angle in the beam
            and frequency [1]
            
        Returns
        -------
        np.array
            TS of a single target at a range of frequencyes [dB re 1m^2]
        """
        TS_m = (
            10 * np.log10(P_rx_e_t_m)
            + 40 * np.log10(r_t)
            + 2 * alpha_m * r_t
            - 10
            * np.log10(
                (p_tx_e * lambda_m**2 * g_theta_t_phi_t_f_t**2) / (16 * np.pi**2)
            )
        )

        return TS_m

    #
    # Chapter IV: VOLUME BACKSCATTERING STRENGTH
    #

    @staticmethod
    def calcSv(
        p_rx_e_n, r_c_n, lambda_f_c, p_tx_e, alpha_f_c, c, tau_eff, psi_f_c, g_0_f_c
    ):
        """
        Calculate volume backscattering strength spectrum at a range step.
        
        Parameters
        ----------
        p_rx_e_n :
            REceived electric power into a matched load [W}]
        r_c_n :
            Range to the centre of sliding window [m]
        lambda_f_c :
            Acoustic wavelength at centre frequency [m]
        p_tx_e : float
            Transmitted electric power [W]
        alpha_f_c : float
            Acoustic absorption at centre frequency [dB/m]
        c : float
           sound speed [m/s] 
        tau_eff : float
            Effective transmit pulse duration [s]
        psi_f_c : float
            Equivalent beam angle at centre frequency [sr]
        g_0_f_c : float
            Transducer gain at centre frequency [1]
            
        Returns
        -------
        np.array
            Sv [dB re 1m^-1]
            
        """
        G = (p_tx_e * lambda_f_c**2 * c * tau_eff * psi_f_c * g_0_f_c**2) / (
            32 * np.pi**2
        )
        Sv_n = (
            10 * np.log10(p_rx_e_n)
            + 20 * np.log10(r_c_n)
            + 2 * alpha_f_c * r_c_n
            - 10 * np.log10(G)
        )

        return Sv_n

    @staticmethod
    def calcPulseCompSphericalSpread(y_pc_n, r_c_n):
        """
        Calculate the spherical spreading compensation.
        
        Parameters
        ----------
        y_pc_n : np.array
            Pulse compressed signal averaged over all transducer sectors [V]
        r_c_n : float
            Range to the centre of of the range volume covered by the sliding 
            window [m]
            
        Returns
        -------
        np.array
            Pulse compressed signal compensated for spherical spreading [Vm]
        """
        y_pc_s_n = y_pc_n * r_c_n

        return y_pc_s_n

    @staticmethod
    def defHanningWindow(c, tau, dr, f_s_dec, N_w = None):
        """
        Calculate the Hann window coefficients.
        
        Length of Hanning window currently chosen as 2^k samples for
        lowest k where 2^k >= 2 * No of samples in pulse
        
        Parameters
        ----------
        c : float
            Sound speed [m/s[]]
        tau : float
            Nominal transmit pulse duration [s]
        dr : float
            Distance between samples [m]
        f_s_dec : float
            Decimated sample rate [Hz]
            
        Returns
        -------
        w_tilde_i : np.array
            Normalised Hann window coefficients [1]
        N_w : float
            Number of samples used in the sliding Hann window [1]
        t_w : float
            Duration of sliding window for calculating volumne backscattering
            strength [s]
        t_w_n : np.array
            Time of each coefficient in `w_tilde_i` [s]    
        """
        if N_w is None:
            L = (c * 2 * tau) / dr  # Number of samples in 2 x pulse duration
            N_w = int(2 ** np.ceil(np.log2(L)))
        t_w_n = np.arange(0, N_w) / f_s_dec
        t_w = N_w / f_s_dec

        w_i = Calculation.hann(N_w)
        w_tilde_i = w_i / (np.linalg.norm(w_i) / np.sqrt(N_w))

        return w_tilde_i, N_w, t_w, t_w_n

    @staticmethod
    def calcDFTforSv(
        y_pc_s_n, w_tilde_i, y_mf_auto_n, N_w, f_m, f_s_dec, r_c_n, step
    ):
        """
        Calculate the normalized DFT of sliding window data.
        
        Parameters
        ----------
        y_pc_s_n :
            Pulse compressed signal compensated for spherical spreading  [Vm]
        w_tilde_i : np.array
            Normalised Hann window coefficients [1]
        y_mf_auto_n : np.array
            Autocorrelation function for the matched filter [1]
        N_w : float
            Number of samples used in the sliding Hann window [1]
        n_f_points : int
            Number of frequencies in `f_m` (not used in this function)
        f_m : np.array
            Frequencies [Hz]
        f_s_dec : float
            Decimated sample rate [Hz]
        r_c_n : float
            Range to the centre of of the range volume covered by the sliding 
            window [m]
        step : float
            Range bin size [m]
        
        Returns
        -------
        Y_pc_v_m_n : np.array
            DFT of the pulse compessed signal from a volume, compensated for
            spreading loss [V]
        Y_mf_auto_m : np.array
            Autocorrection function for the matched filter [1]
        Y_tilde_pc_v_m_n : np.array
            DFT of the pulse compressed signal from a volume normalised by the DFT
            of the reduced autocorrelation function for the matched filter,
            compensated for spreading loss [1]
        svf_range : 
            Ranges at which the returned variables are at [m]
        """
        # Prepare for append
        Y_pc_v_m_n = []
        Y_tilde_pc_v_m_n = []
        svf_range = []
    
        # DFT of auto correlation function of the matched filter signal
        _Y_mf_auto_m = np.fft.fft(y_mf_auto_n, n=N_w)
        Y_mf_auto_m = Calculation.freqtransf(_Y_mf_auto_m, f_s_dec, f_m)

        min_sample = 0  # int(r0 / dr)
        max_sample = len(y_pc_s_n)  # int(r1 / dr)

        bin_start_sample = min_sample
        bin_stop_sample = bin_start_sample + N_w

        while bin_stop_sample < max_sample:
            # Windowed data
            yspread_bin = w_tilde_i * y_pc_s_n[bin_start_sample:bin_stop_sample]

            # Range for bin
            bin_center_sample = int((bin_stop_sample + bin_start_sample) / 2)
            bin_center_range = r_c_n[bin_center_sample]
            svf_range.append(bin_center_range)

            # DFT of windowed data
            _Y_pc_v_m = np.fft.fft(yspread_bin, n=N_w)
            Y_pc_v_m = Calculation.freqtransf(_Y_pc_v_m, f_s_dec, f_m)

            # Normalized DFT of windowed data
            Y_tilde_pc_v_m = Y_pc_v_m / Y_mf_auto_m

            # Append data
            Y_pc_v_m_n.append([Y_pc_v_m])
            Y_tilde_pc_v_m_n.append([Y_tilde_pc_v_m])

            # Next range bin
            bin_start_sample += step
            bin_stop_sample = bin_start_sample + N_w

        svf_range = np.array(svf_range)

        return Y_pc_v_m_n, Y_mf_auto_m, Y_tilde_pc_v_m_n, svf_range

    @staticmethod
    def freqtransf(FFTvecin, fsdec, fvec=None):
        """
        Shift FFT frequencies for specified frequencies.
        
        Parameters
        ----------
        FFTvecin : np.array
            FFT data from decimated frequencies
        fsdec : float
            Decimated sampling frequency [Hz]
        fvec : np.array
            Specified frequencies [Hz]
            
            From calibration data. ("CAL": {"frequencies")
                      If no calibration generate freq vector starting from 
                       f0 to f1 with same number of points as in calibration data
                    
        Returns
        -------
        float
            Vector with corrected frequencies [Hz]
        """

        nfft = len(FFTvecin)
        idxtmp = np.floor(fvec / fsdec * nfft).astype("int")
        idx = np.mod(idxtmp, nfft)

        return FFTvecin[idx]

    @staticmethod
    def calcPowerFreqSv(Y_tilde_pc_v_m_n, N_u, z_rx_e, z_td_e):
        """
        Calculate the received power spectrum for the sliding window.
        
        Parameters
        ----------
        Y_tilde_pc_v_m_n : np.array
            DFT of the pulse compressed signal from a volume normalised by the DFT
            of the reduced autocorrelation function for the matched filter,
            compensated for spreading loss [1]
        N_u : int
            Number of transducer sectors/receiver channels
        z_td_e : float
            Transducer sector electric impedance [Ω]
        z_rx_e : float
            Receiver electric impedance [Ω]
        
        Returns
        -------
        np.array
            DFT of the received electric power in a matched load for the signal
            from a volume [Wm^2]
        """
        
        # Initialize list of power values by range
        P_rx_e_v_m_n = []

        # Impedances
        Z = (np.abs(z_rx_e + z_td_e) / np.abs(z_rx_e)) ** 2 / np.abs(z_td_e)

        # Loop over list of FFTs along range
        for Y_tilde_pc_v_m in Y_tilde_pc_v_m_n:
            P_rx_e_v_m = N_u * (np.abs(Y_tilde_pc_v_m) / (2 * np.sqrt(2))) ** 2 * Z
            # Append power to list
            P_rx_e_v_m_n.append(P_rx_e_v_m)

        return P_rx_e_v_m_n

    @staticmethod
    def calcSvf(P_rx_e_t_m_n, alpha_m, p_tx_e, lambda_m, t_w,
                psi_m, g_0_m, c, svf_range):
        """
        Calculate Sv as a function of frequency.
        
        Parameters
        ----------
        P_rx_e_t_m_n : np.array
            DFT of the received electric power [W]
        alpha_m : float
            Acoustic absorption [dB/m]
        p_tx_e : float
            Transmitted electric power [W]
        lambda_m : float
            Acoustic wavelength [m]
        t_w : float
            Sliding window duration [s]
        psi_m : float
            Equivalent beam angle [sr]
        g_0_m : float
            Transducer gain [dB]
        c : float
            Speed of sound [m/s]
        svf_range : np.array
            Range [m]
            
        Returns
        -------
        np.array
            Sv(f) [dB re 1 m^-1]
        """
        
        # Initialize list of Svf by range
        Sv_m_n = np.empty([len(svf_range), len(alpha_m)], dtype=float)

        G = (p_tx_e * lambda_m**2 * c * t_w * psi_m * g_0_m**2) / (32 * np.pi**2)
        n = 0

        # Loop over list of power values along range
        for P_rx_e_t_m in P_rx_e_t_m_n:
            Sv_m = (
                10 * np.log10(P_rx_e_t_m)
                + 2 * alpha_m * svf_range[n]
                - 10 * np.log10(G)
            )

            # Add to array
            Sv_m_n[n, ] = Sv_m
            n += 1

        return Sv_m_n

    @staticmethod
    def calcpsi(psi_f_n, f_n, f_m):
        """
        Calculate psi at given frequency.
        
        Parameters
        ----------
        psi_f_n : float
            Psi at nominal frequency [sr]
        f_n : float
            Nominal frequency [Hz]
        f_m : float
            Frequency to calculate psi at [Hz]
            
        Returns
        -------
        float 
            Psi at frequency `f_m` [sr]
        """
        
        return psi_f_n * (f_n / f_m) ** 2
    

    @staticmethod
    def calcalpha(data, f_m):

        def alphaFG(c, pH, T, D, S,f): # requires sound speed (m/s), pH, temp(C), depth(m), salinity(ppt), and nominal frequency(kHz)
            # Attenuation Coefficient is based on Francois and Garrison, 1982 - "Sound absorption based on ocean measurements.
            # Boric Acid Contribution, P1 = 1. This is buried in echolab.simrad_calibration and I can't figure out how to apply to the 
            # full range of frequencies and not the nominal frequency parameter in the raw data
            A1=((8.86/c)*(10**(0.78*pH-5)))
            f1=((2.8*((S/35)**0.5))*(10**(4-(1245/(T+273)))))
            # MgSO4 Contribution
            A2=((21.44*(S/c))*(1+(0.025*T)))
            P2=(1-(1.37*(10**-4)*D)+(6.2*(10**-9)*(D**2)))
            f2=((8.17*(10**(8-(1990/(T+273)))))/(1+.0018*(S-35)))
            # Pure water contribution, where A3 is temperature dependent
            if T > 20:
                A3=((3.964*(10**-4))-(1.146*(10**-5)*T)+(1.45*(10**-7)*(T**2))-(6.5*(10**-10)*(T**3)))
            else:
                A3=((4.937*(10**-4))-(2.59*(10**-5)*T)+(9.11*(10**-7)*(T**2))-(1.5*(10**-8)*(T**3)))
            P3=((1-(3.83*(10**-5)*D)) + (4.9*(10**-10)*(D**2)))
            # Calculate and return Alpha
            alpha = (((f**2)*A1*f1)/(((f1**2)) + (f**2)))+ ((A2*P2*f2*(f**2))/((f2**2) + (f**2))) + (A3*P3*(f**2))
            return alpha
        
        env = data.environment[data.environment != np.array(None)][0]
        alpha_m = np.array([alphaFG(env['sound_speed'],env['acidity'],env['temperature'],env['depth'],env['salinity'],nomf/1000)/1000 for nomf in f_m])
        
        return alpha_m