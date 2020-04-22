import numpy as np
import cv2
import matplotlib.pyplot as plt


class Time:
    def __init__(self, wave, dt):
        self.wave = wave
        self.n = len(self.wave)
        self.mean = np.mean(self.wave)
        self.std = np.std(self.wave)
        self.dt = dt
        self.df = 1 / self.n / self.dt

    def set_wave(self, wave, dt):
        self.wave = wave
        self.n = len(self.wave)
        self.mean = np.mean(self.wave)
        self.std = np.std(self.wave)
        self.dt = dt
        self.df = 1 / self.n / self.dt

    def set_dt(self, dt):
        self.dt = dt
        self.df = 1 / self.n / self.dt

    @property
    def fft(self):
        return np.fft.fft(self.wave)

    @property
    def fft_abs(self):
         return np.abs(self.fft) / self.n * 2

    def image(self, width):
        image = np.tile(self.wave.astype(np.uint8), (width, 1))
        return image

    @property
    def t_sq(self):
        return np.linspace(0, self.dt * self.n, self.n)

    @property
    def fq_sq(self):
        return np.linspace(0, 1.0 / self.dt, self.n)

    def ss_curve(self, ps):
        return 5.05 * (np.exp(-0.138 * (1 / 180 * np.pi) * 600 * self.fq_sq / ps) * (
                    1 - np.exp(-0.1 * (1 / 180 * np.pi) * 600 * self.fq_sq / ps)))


class Freq:
    def __init__(self, fft, df):
        self.fft = fft
        self.n = len(self.fft)
        self.df = df
        self.dt = 1 / self.n / self.df

    def set_f(self, fft):
        self.fft = fft
        self.n = len(self.fft)
        self.dt = 1 / self.n / self.df

    def set_df(self, df):
        self.n = len(self.fft)
        self.df = df
        self.dt = 1 / self.n / self.df

    @property
    def fq_sq(self):
        return np.linspace(0, 1.0 / self.dt, self.n)

    @property
    def fft_abs(self):
        return np.abs(self.fft) / self.n * 2

    @property
    def ifft(self):
        return np.fft.ifft(self.fft)

    def inv_wave(self, *args):
        return self.ifft.real
