import cv2
import matplotlib.patches as patches
import matplotlib.pylab as plt
from matplotlib.widgets import RadioButtons , RectangleSelector, Slider
import numpy as np
import imageanalyzer2
import math
import copy
# import PySimpleGUI as sg


def division_quotient(max_num):
    """
    グラフのx軸の間隔をちょうど良い感じにする目的の関数。
    2,5,10の倍数のいずれかの単位でmax_numを５〜１０分割する
    戻りは分割単位
    """
    delta = 0
    mod = max_num
    for i in (5 , 6 , 7 , 8):  # グラフをざっくり分割する個数
        x = max_num / i
        k = int ( math.log10 ( x ) )
        if x < 1: k -= 1
        for j in (1 , 2 , 5):  # 2,5,10の倍数で等分する
            tmpdelta = j * 10 ** k
            tmpmod = max_num - tmpdelta * i
            if abs ( mod ) > abs ( tmpmod ):  # あまりが小さいのものを選ぶ
                mod = tmpmod
                delta = tmpdelta
    return delta


class Simulation ():
    def __init__(self , file , dt , ps=200, isimage=True):
        self.file = file
        self.dt = dt
        self.ps = ps
        self.datanum = None
        self.f_range = 300
        self.amp_max = 5
        self.wave_cont = 1
        self.isimage = isimage
        self.ref_timeobj = imageanalyzer2.Time ()
        self.sim_timeobj = imageanalyzer2.Time ()
        self.ref_fqobj = imageanalyzer2.Freq ()
        self.sim_fqobj = imageanalyzer2.Freq ()
        self.ref_wave = None
        self.t_sq = self.ref_timeobj.t_sq
        self.fq_sq = self.ref_fqobj.fq_sq
        self.x_a , self.y_a , self.x_b , self.y_b = None , None , None , None
        self.patch = patches.Rectangle ( xy=(0,0), width=1 , height=1, fill=False )
        self.labels = ('pixel','mm','sec')
        self.label = self.labels[0]
        self.filter = None
        self.selected_area = None  # or 'fft'
        self.set_fig ()
        self.set_ax ()
        self.set_frame ()
        self.set_widgets()

    def set_fig(self):
        self.fig , (self.ax_original_data , self.ax_ref_image , self.ax_sim_image , \
                    self.ax_wave , self.ax_fft) = plt.subplots ( nrows=5 , figsize=(8 , 8) )
        self.fig.suptitle('Banding Designer', fontsize=16)

    def set_ax(self):
        self.line_refwave , = self.ax_wave.plot ( [] , [] , 'green' , \
                                                  linestyle="dashed" , linewidth=0.5 , label='ref' )
        self.line_reffft , = self.ax_fft.plot ( [] , [] , 'green' , \
                                                linestyle="dashed" , linewidth=0.5 , label='ref' )
        self.line_simwave , = self.ax_wave.plot ( [] , [] , 'black' , \
                                                  linewidth=1 , label='sim' )
        self.line_simfft , = self.ax_fft.plot ( [] , [] , 'black' , \
                                                linewidth=1 , label='sim' )
        self.line_vtffft , = self.ax_fft.plot ( [] , [] , 'gray' , \
                                                linestyle="dashed" , linewidth=0.5 , label='vtf' )

    def set_frame(self):
        self.ax_original_data.set_title ( 'original image' )
        self.ax_ref_image.set_title ( 'ref image' )
        self.ax_sim_image.set_title ( 'sim image' )
        self.ax_ref_image.set_yticks([])
        self.ax_sim_image.set_yticks([])
        self.ax_ref_image.set_xticks([])
        self.ax_sim_image.set_xticks([])
        for spine in ('left' , 'right' , 'top' , 'bottom'):
            self.ax_ref_image.spines[spine].set_color ( "g" )
            self.ax_ref_image.spines[spine].set_linewidth ( 1 )
            self.ax_ref_image.spines[spine].set_linestyle ( "dashed" )
            self.ax_sim_image.spines[spine].set_color ( "black" )
            self.ax_sim_image.spines[spine].set_linewidth ( 2 )

        self.ax_wave.set_xlabel ( 'sec' )
        self.ax_wave.set_ylabel ( 'L*' )
        self.ax_wave.grid ( True )
        self.ax_wave.legend (loc='upper right')
        self.ax_fft.set_xlabel ( 'Hz' )
        self.ax_fft.set_ylabel ( 'amp' )
        self.ax_fft.grid ( True )
        self.ax_fft.legend (loc='upper right')

    def set_x_scale(self, label):
        self.label = label
        if self.datanum is not None:
            # ref/sim_image xtick設定
            xticks , xlabels  = self.set_image_x_array(self.datanum )
            self.set_image_x_sclae ( self.ax_ref_image , xticks , xlabels )
            self.set_image_x_sclae ( self.ax_sim_image , xticks , xlabels )
        # original_image xtick設定
        xticks , xlabels  = self.set_image_x_array(self.original_datanum )
        self.set_image_x_sclae ( self.ax_original_data , xticks , xlabels )
        self.fig.canvas.draw()

    def set_image_x_array(self, data_num):
        max_data_num, delta_pix, dq = None, None, None
        if self.label == 'pixel':
            max_data_num = data_num
            dq = division_quotient ( max_data_num )
            delta_pix = dq
        elif self.label == 'sec':
            max_data_num = data_num * self.dt
            dq = division_quotient ( max_data_num )
            delta_pix = int( dq / self.dt )
        elif self.label == 'mm':
            max_data_num = data_num * self.dt * self.ps
            dq = division_quotient ( max_data_num )
            delta_pix = int ( dq / self.dt / self.ps )
        ticks = np.arange(0, data_num, delta_pix)
        labels  = np.round( np.arange( 0 , max_data_num , dq ) , 1 )
        return ticks, labels

    def set_image_x_sclae(self, ax, xticks, xlabels):
        ax.set_xticks (xticks)
        ax.set_xticklabels (xlabels)

    def set_widgets(self):
        axcolor = 'lightgoldenrodyellow'
        ra_xs = plt.axes([0.8, 0.9, 0.1, 0.1], facecolor=axcolor)
        self.radio_x_scale = RadioButtons(ra_xs, self.labels, activecolor='green' )

        ax_freq = plt.axes ( [0.2 , 0.02 , 0.6 , 0.02] , facecolor=axcolor )
        self.s_freq = Slider ( ax_freq , 'max freq' , 50 , 2000 , valinit=self.f_range \
                               , valstep=50 , valfmt='%1.0f', color='green' )

        ax_amp = plt.axes ( [0.02 , 0.12 , 0.02 , 0.12] , facecolor=axcolor )    # axes([左, 下, 幅, 高さ])
        self.s_amp = Slider ( ax_amp , 'amp' , 1 , 10, valinit=self.amp_max \
                 ,orientation='vertical', valstep=1, valfmt='%1.0f', color='green' )

        ax_cont = plt.axes ( [0.02 , 0.28 , 0.02 , 0.12] , facecolor=axcolor )    # axes([左, 下, 幅, 高さ])
        self.s_cont = Slider ( ax_cont , 'cont' , 0.2 , 1.8, valinit=self.wave_cont \
                 ,orientation='vertical', valstep=0.1, valfmt='%1.1f', color='green' )

    def onselect(self , eclick , erelease):
        pass  # 座標を返すこともできるが、on,offの座標が逆転するバグがあるので使わない

    def onclick(self , event):
        if event.inaxes == self.ax_original_data:
            if event.button == 1:
                self.x_a , self.y_a = event.xdata , event.ydata
                self.selected_area = 'origin'
            elif event.button == 2 or event.button == 3:
                self.patch.set_visible ( not self.patch.get_visible () )
                self.fig.canvas.draw ()

        elif event.inaxes == self.ax_fft:
            if event.button == 1:
                self.x_a , self.y_a = event.xdata , event.ydata
                self.selected_area = 'fft'
            elif event.button == 2:  # reset
                self.sim_timeobj.set_wavedt ( self.ref_wave , self.dt )
                self.sim_fqobj.set_f ( self.ref_timeobj.fft )
                self.sim_draw ()
                self.fig.canvas.draw ()
            elif event.button == 3:  # undo
                pre_filter = 1 / self.filter
                self.sim_fqobj.set_f(pre_filter * self.sim_fqobj.fft)
                self.sim_timeobj.set_wavedt ( self.sim_fqobj.inv_wave () , self.dt )
                self.filter = pre_filter
                self.sim_draw ()
                self.fig.canvas.draw ()

    def offclick(self, event):
        self.x_b , self.y_b = event.xdata , event.ydata
        if event.button == 1:
            if self.selected_area == 'origin' and event.inaxes == self.ax_original_data:
                self.ref_setting ()
            elif self.selected_area == 'fft' and event.inaxes == self.ax_fft:
                self.fft_change ()
            self.selected_area = None
            self.fig.canvas.draw ()

    def update(self, val):
        self.f_range = self.s_freq.val
        self.amp_max = self.s_amp.val
        self.ax_fft.set_xlim ( 0 , self.s_freq.val )
        self.ax_fft.set_ylim(0, self.s_amp.val)
        self.fig.canvas.draw_idle ()

    def update_cont(self, val):
        self.wave_cont = self.s_cont.val
        prop_fft = self.sim_fqobj.fft /self.ref_fqobj.fft
        self.ref_timeobj.set_wave (
            (self.ref_wave - np.mean(self.ref_wave) ) * self.wave_cont + np.mean(self.ref_wave) )
        self.ref_fqobj.set_f( self.ref_timeobj.fft)
        self.sim_fqobj.set_f(prop_fft * self.ref_fqobj.fft)  # 記憶していたFFT変更比率を掛ける
        self.sim_timeobj.set_wave ( self.sim_fqobj.inv_wave () )
        self.ref_draw()
        self.sim_draw()
        self.fig.canvas.draw_idle ()

    def ref_setting(self):
        self.patch.set_visible(False)
        if self.isimage:    self.set_original_image()
        else:   self.set_original_wave()
        self.trim_original_data()
        self.set_initial_wave ()
        self.ref_draw ()
        self.sim_draw ()
        self.set_initial_axis ()

    def trim_original_data(self):
        top = int ( min ( self.x_a , self.x_b ) )
        end = int ( max ( self.x_a , self.x_b ) )
        left = int ( min ( self.y_a , self.y_b ) )
        right = int ( max ( self.y_a , self.y_b ) )
        self.patch = patches.Rectangle ( xy=(top , left) , width=(end - top) , height=(right - left) , fill=False )
        self.ax_original_data.add_patch ( self.patch )
        if self.isimage:
            trimmed_image = self.raw_image[top: end , left: right]
            self.ref_wave = np.mean ( trimmed_image , axis=1 )
        else:
            self.ref_wave = self.original_wave[top: end]

    def set_initial_wave(self):
        self.ref_timeobj.set_wavedt (
            (self.ref_wave - np.mean(self.ref_wave) ) * self.wave_cont + np.mean(self.ref_wave), self.dt )
        self.datanum = self.ref_timeobj.n
        self.ref_fqobj.set_fdf ( self.ref_timeobj.fft , self.ref_timeobj.df )
        self.sim_fqobj = copy.copy(self.ref_fqobj)
        self.sim_timeobj.set_wavedt ( self.sim_fqobj.inv_wave () , self.dt )
        self.t_sq = self.ref_timeobj.t_sq
        self.fq_sq = self.ref_fqobj.fq_sq

    def set_initial_axis(self):
        self.ax_wave.set_xlim ( 0 , self.datanum * self.ref_timeobj.dt )
        self.ax_wave.set_ylim (min(self.ref_wave)*0.95, max(self.ref_wave)*1.05)
        self.ax_fft.set_xlim ( 0 , self.f_range )
        self.ax_fft.set_ylim ( 0 , self.amp_max )
        self.set_x_scale(self.label)

    def ref_draw(self):
        self.line_refwave.set_data ( self.t_sq , self.ref_timeobj.wave )
        self.line_reffft.set_data ( self.fq_sq , self.ref_fqobj.fft_abs )
        self.ax_ref_image.imshow ( self.ref_timeobj.image ( 1 ) , 'gray' , vmin=0 , vmax=255 , aspect='auto' )
        self.line_vtffft.set_data ( self.fq_sq , self.ref_timeobj.ss_curve(self.ps))

    def sim_draw(self):
        self.line_simfft.set_data ( self.fq_sq , self.sim_fqobj.fft_abs )
        self.line_simwave.set_data ( self.t_sq , self.sim_timeobj.wave )
        self.ax_sim_image.imshow ( self.sim_timeobj.image ( 1 ) , 'gray' , vmin=0 , vmax=255 , aspect='auto' )

    def fft_change(self):
        self.pre_prop_fft = self.sim_fqobj.fft / self.ref_fqobj.fft
        self.set_fft_filter ()
        self.sim_fqobj.set_f ( self.sim_fqobj.fft * self.filter )
        self.sim_timeobj.set_wavedt ( self.sim_fqobj.inv_wave () , self.dt )
        self.line_simwave.set_data ( self.t_sq , self.sim_timeobj.wave )
        self.line_simfft.set_data ( self.fq_sq , self.sim_fqobj.fft_abs )
        self.ax_sim_image.imshow ( self.sim_timeobj.image ( 1 ) , 'gray' , vmin=0 , vmax=255 , aspect='auto' )

    def set_fft_filter(self):
        self.filter = np.ones ( self.datanum)
        prop = self.y_b / self.y_a
        out_freq = max ( self.x_a , self.x_b )
        in_freq = min ( self.x_a , self.x_b )
        if in_freq > out_freq:
            tmp = in_freq
            in_freq = out_freq
            out_freq = tmp
        self.filter[(self.fq_sq > in_freq) & (self.fq_sq < out_freq)] = prop
        self.filter[(self.fq_sq > self.fq_sq[-1] - out_freq) & \
                    (self.fq_sq < self.fq_sq[-1] - in_freq)] = prop

    def set_original_image(self):
        self.raw_image = cv2.imdecode ( np.fromfile ( self.file , np.uint8 ) , 0 )
        self.original_datanum = len(self.raw_image)
        self.ax_original_data.imshow ( self.raw_image.T , 'gray' , vmin=0 , vmax=255 , aspect='auto' )

    def set_original_wave(self):
        self.original_wave = np.loadtxt ( 'sinwave.txt' , dtype='uint8' )
        self.org_timeobj = imageanalyzer2.Time ( wave=self.original_wave , dt=self.dt )
        self.original_datanum = self.org_timeobj.n
        dot_sq = np.arange(0, self.original_datanum )
        self.line_origilal_wave , = self.ax_original_data.plot (dot_sq, self.org_timeobj.wave, 'black' , \
                                                                 linewidth=1 , label='origin' )
        self.ax_original_data.grid ( True )
        self.ax_original_data.set_ylim (min(self.org_timeobj.wave)*0.9, max(self.org_timeobj.wave)*1.1)
        self.ax_original_data.set_ylabel ( 'L*' )

    def run(self):
        if self.isimage:
            self.ax_original_data.set_title ( 'original image' )
            self.set_original_image ()
        else:
            self.ax_original_data.set_title ( 'original wave' )
            self.set_original_wave ()

        self.radio_x_scale.on_clicked(self.set_x_scale)
        RS_im = RectangleSelector ( self.ax_original_data , self.onselect , drawtype='box' , useblit=True )
        RS_fft = RectangleSelector ( self.ax_fft , self.onselect , drawtype='box' , useblit=True )
        self.fig.canvas.mpl_connect ( 'button_press_event', self.onclick )
        self.fig.canvas.mpl_connect ( 'button_release_event', self.offclick )
        self.s_freq.on_changed( self.update )
        self.s_amp.on_changed( self.update )
        self.s_cont.on_changed( self.update_cont )

        plt.show ()


if __name__ == '__main__':
    dt = 0.0005
    # file = 'sinwave.txt'
    file = 'ref.tif'
    sim = Simulation ( file, dt , ps=200, isimage=1)
    sim.run ()
