import cv2
import matplotlib.patches as patches
import matplotlib.pylab as plt
from matplotlib.widgets import RadioButtons, RectangleSelector
import numpy as np
import imageanalyzer2
import PySimpleGUI as sg

class Simuration():
    def __init__(self, file, dt, isimage=True):
        self.file = file
        self.dt = dt
        # self.ps = 226
        self.f_range = 100
        self.amp_max = 5
        self.isimage = isimage
        self.ref_timeobj = imageanalyzer2.Time(np.array([0, 1]), dt)
        self.sim_timeobj = imageanalyzer2.Time(np.array([0, 1]), dt)
        self.ref_fqobj = imageanalyzer2.Freq(np.array([0, 1]), 1)
        self.sim_fqobj = imageanalyzer2.Freq(np.array([0, 1]), 1)
        self.pre_fft = None
        self.ref_wave = None
        self.t_sq = self.ref_timeobj.t_sq
        self.fq_sq = self.ref_fqobj.fq_sq
        self.x_a, self.y_a, self.x_b, self.y_b = None, None, None, None
        self.filter = None
        self.area_image = False
        self.area_fft = False
        self.set_fig()
        self.set_ax()
        self.set_frame()

    def set_fig(self):
        self.fig, (self.ax_original_image, self.ax_ref_image, self.ax_sim_image,\
                   self.ax_wave, self.ax_fft) = \
        plt.subplots(nrows=5, figsize=(8,10))

    def set_ax(self):
        self.line_refwave, = self.ax_wave.plot([], [], 'green',\
                                               linestyle="dashed", linewidth=0.5, label='ref')
        self.line_reffft, =  self.ax_fft.plot([], [], 'green',\
                                               linestyle="dashed", linewidth=0.5, label='ref')
        self.line_simwave, = self.ax_wave.plot([], [], 'black',\
                                               linewidth=1, label='SIM')
        self.line_simfft, =  self.ax_fft.plot([], [], 'black',\
                                               linewidth=1, label='SIM')

    def ref_draw(self):
        self.line_refwave.set_data(self.x_data, self.ref_timeobj.wave)
        self.line_reffft.set_data(self.fq_sq, self.ref_fqobj.fft_abs)
        self.ax_ref_image.imshow(self.ref_timeobj.image(1), 'gray', vmin=0, vmax=255, aspect='auto')

    def sim_draw(self):
        self.line_simfft.set_data(self.fq_sq, self.sim_fqobj.fft_abs)
        # self.line_simwave.set_data(self.x_data, self.sim_timeobj.wave)
        self.ax_sim_image.imshow(self.sim_timeobj.image(1), 'gray', vmin=0, vmax=255, aspect='auto')

    def set_frame(self):
        self.ax_original_image.set_title('original image')
        self.ax_ref_image.set_title('ref image')
        # todo refactoring
        self.ax_ref_image.spines["right"].set_color("g")
        self.ax_sim_image.set_title('sim image')
        self.ax_wave.grid(True)
        self.ax_fft.grid(True)
        self.ax_wave.set_xlabel('sec')
        self.ax_wave.set_ylabel('L*')
        self.ax_fft.set_xlabel('Hz')
        self.ax_fft.set_ylabel('amp')
        self.ax_wave.legend()
        self.ax_fft.legend()

    def fft_filter(self, in_freq, out_freq, prop):
        self.filter = np.ones(self.ref_timeobj.n)
        if in_freq > out_freq:
            tmp = in_freq
            in_freq = out_freq
            out_freq = tmp
        # todo refactoring
        self.filter[(self.fq_sq > in_freq) & (self.fq_sq < out_freq)] = prop
        self.filter[(self.fq_sq > self.fq_sq[-1] - out_freq) & \
                    (self.fq_sq < self.fq_sq[-1] - in_freq)] = prop

    def onselect(self, eclick, erelease):
        pass

    def onclick(self, event):
        if event.inaxes == self.ax_original_image:
            if event.button == 1:
                self.x_a, self.y_a = event.xdata, event.ydata
                self.area_image = True
            elif event.button == 2 or event.button == 3:
                self.reselect_area()
                self.fig.canvas.draw()

        elif event.inaxes == self.ax_fft:
            if event.button == 1:
                self.x_a, self.y_a = event.xdata, event.ydata
                self.area_fft = True
            elif event.button == 2:     # reset
                self.sim_timeobj.set_wave(self.ref_wave)
                self.sim_fqobj.set_f(self.ref_timeobj.fft)
                self.sim_draw()
                self.fig.canvas.draw()

            elif event.button == 3:     # undo
                tmp_fft = self.sim_fqobj.fft
                self.sim_fqobj.set_f(self.pre_fft)
                self.pre_fft = tmp_fft
                self.sim_timeobj.set_wave(self.sim_fqobj.inv_wave())
                self.sim_draw()
                self.fig.canvas.draw()

    def reselect_area(self):
        self.ax_original_image.cla()
        self.set_original_image()
        self.ax_ref_image.cla()
        self.ax_sim_image.cla()
        self.line_refwave.remove()
        self.line_reffft.remove()
        self.line_simwave.remove()
        self.line_simfft.remove()
        self.set_frame()

    def offclick(self, event):
        self.x_b, self.y_b = event.xdata, event.ydata
        if event.button == 1:
            if self.area_image and event.inaxes == self.ax_original_image:
                self.ref_setting()
            elif self.area_fft and event.inaxes == self.ax_fft:
                self.fft_change()
            self.area_image = False
            self.area_fft = False
            self.fig.canvas.draw()

    def ref_setting(self):
        self.ax_original_image.cla()
        self.set_original_image()
        self.trim_ref_image()
        self.ref_wave = np.mean(self.trimmed_image, axis=1)
        self.set_initial_wave()
        self.ref_draw()
        self.sim_draw()
        self.set_initial_axis()

    def trim_ref_image(self):
        top = int(min(self.x_a, self.x_b))
        end = int(max(self.x_a, self.x_b))
        left = int(min(self.y_a, self.y_b))
        right = int(max(self.y_a, self.y_b))
        r = patches.Rectangle(xy=(top, left), width=(end-top), height=(right-left), fill=False)
        self.ax_original_image.add_patch(r)
        self.trimmed_image = self.raw_image[top: end, left: right]

    def set_initial_wave(self):
        self.ref_timeobj.set_wave(self.ref_wave)
        self.ref_fqobj.set_f(self.ref_timeobj.fft)
        self.ref_fqobj.set_df(self.ref_timeobj.df)
        self.sim_fqobj = self.ref_fqobj
        self.sim_timeobj.set_wave(self.sim_fqobj.inv_wave())
        self.t_sq = self.ref_timeobj.t_sq
        self.x_data = self.t_sq
        self.fq_sq = self.ref_fqobj.fq_sq

    def set_initial_axis(self):
        self.ax_wave.set_xlim(0, self.ref_timeobj.n * self.ref_timeobj.dt)
        self.ax_wave.set_ylim(0, 255)
        self.ax_fft.set_xlim(0, self.f_range)
        self.ax_fft.set_ylim(0, self.amp_max)

    def fft_change(self):
        prop = self.y_b / self.y_a
        out_freq = max(self.x_a, self.x_b)
        in_freq = min(self.x_a, self.x_b)
        self.pre_fft = self.sim_fqobj.fft
        self.fft_filter(in_freq, out_freq, prop)
        self.sim_fqobj.set_f(self.pre_fft * self.filter)
        self.sim_timeobj.set_wave(self.sim_fqobj.inv_wave())
        self.line_simwave.set_data(self.x_data, self.sim_timeobj.wave)
        self.line_simfft.set_data(self.fq_sq, self.sim_fqobj.fft_abs)
        self.ax_sim_image.imshow(self.sim_timeobj.image(1),'gray', vmin=0, vmax=255, aspect='auto')

    def set_original_image(self):
        self.raw_image = cv2.imdecode(np.fromfile(self.file, np.uint8), 0)
        self.ax_original_image.imshow(self.raw_image.T, 'gray', vmin=0, vmax=255, aspect='auto')

    def run(self):
        if self.isimage:
            self.set_original_image()

        else: pass

        RS_im = RectangleSelector(self.ax_original_image, self.onselect, drawtype='box', useblit=True)
        RS_fft = RectangleSelector(self.ax_fft,  self.onselect, drawtype='box', useblit=True)
        self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.fig.canvas.mpl_connect('button_release_event', self.offclick)

        plt.show()

if __name__ =='__main__':
    dt = 0.0005
    sim = Simuration(('ref.tif'), dt, isimage=1)
    sim.run()