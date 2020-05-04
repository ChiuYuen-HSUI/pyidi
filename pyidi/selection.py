import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from multiprocessing import Pool
from tqdm import tqdm

from . import pyidi

class ManualROI:
    """Manual ROI selection."""

    def __init__(self, video, roi_size, single=False, verbose=0):
        """Manually select region of interest.
        
        :param video: parent object
        :type video: object
        :param roi_size: size of region of interest (dy, dx)
        :type roi_size: tuple, list
        :param single: if True, ONLY ONE ROI can be selected, defaults to False
        :type single: bool, optional
        :param verbose: Show text, defaults to 0
        :type verbose: int, optional
        """
        self.roi_size = roi_size
        self.image = video.mraw[0]
        self.verbose = verbose

        # Tkinter root and matplotlib figure
        root = tk.Tk()
        root.title('Pick point')
        fig = Figure(figsize=(15, 7))
        ax = fig.add_subplot(111)
        ax.grid(False)
        ax.imshow(self.image, cmap='gray')
        plt.show()

        # Embed figure in tkinter winodw
        canvas = FigureCanvasTkAgg(fig, root)
        canvas.get_tk_widget().pack(side='top', fill='both', expand=1)
        NavigationToolbar2Tk(canvas, root)

        if self.verbose:
            print('SHIFT + LEFT mouse button to pick a pole.\nRIGHT mouse button to erase the last pick.')


        self.point = [[], []]
        line, = ax.plot(self.point[1], self.point[0], 'r.')
        
        self.rectangles = []
        self.rectangles.append(patches.Rectangle((0, 0), 10, 10, fill=False, alpha=0))
        ax.add_patch(self.rectangles[-1])        
        
        self.shift_is_held = False
        def on_key_press(event):
            """Function triggered on key press (shift)."""
            if event.key == 'shift':
                self.shift_is_held = True

        def on_key_release(event):
            """Function triggered on key release (shift)."""
            if event.key == 'shift':
                self.shift_is_held = False

        def onclick(event):
            if event.button == 1 and self.shift_is_held:
                if event.xdata is not None and event.ydata is not None:
                    if single:
                        self.point[0] = [int(np.round(event.ydata))]
                        self.point[1] = [int(np.round(event.xdata))]
                    else:
                        self.point[0].append(int(np.round(event.ydata)))
                        self.point[1].append(int(np.round(event.xdata)))
                    if self.verbose:
                        print(f'y: {np.round(event.ydata):5.0f}, x: {np.round(event.xdata):5.0f}')

            elif event.button == 3 and self.shift_is_held and not single:
                if self.verbose:
                    print('Deleted the last point...')
                del self.point[1][-1]
                del self.point[0][-1]
                del self.rectangles[-1]

            line.set_xdata(self.point[1])
            line.set_ydata(self.point[0])
            
            if self.point[0]:
                [p.remove() for p in reversed(ax.patches)]
                self.rectangles = []
                for i, (p0, p1) in enumerate(zip(self.point[0], self.point[1])):
                    self.rectangles.append(patches.Rectangle((p1 - self.roi_size[1]//2, p0 - self.roi_size[0]//2), 
                                                    self.roi_size[1], self.roi_size[0], fill=False, color='r', linewidth=2))
                    ax.add_patch(self.rectangles[-1])

            fig.canvas.draw()

        def handle_close(event):
            """On closing."""
            self.points = np.asarray(self.point).T
            if self.verbose:
                for i, point in enumerate(self.polygon):
                    print(f'{i+1}. point: x ={point[1]:5.0f}, y ={point[0]:5.0f}')

        # Connecting functions to event manager
        fig.canvas.mpl_connect('key_press_event', on_key_press)
        fig.canvas.mpl_connect('key_release_event', on_key_release)
        fig.canvas.mpl_connect('button_press_event', onclick)
        # on closing the figure
        fig.canvas.mpl_connect('close_event', handle_close)

        root.mainloop()
    

class GridOfROI:
    """
    Automatic simple ROI grid generation.

    Different from RegularROIGrid in that it gets a regular grid and only
    then checks if all points are inside polygon. This yields a more regular
    and full grid. Does not contain sssig filter.
    """
    def __init__(self, video=None, roi_size=(7, 7), noverlap=0, verbose=0):
        """
        
        :param video: parent object of video
        :type video: object
        :param roi_size: Size of the region of interest (y, x), defaults to (7, 7)
        :type roi_size: tuple, list, optional
        :param noverlap: number of pixels that overlap between neighbouring ROIs
        :type noverlap: int, optional
        :param sssig_filter: minimum value of SSSIG that the roi must have, defaults to None
        :type sssig_filter: None, float, optional
        :param verbose: Show text, defaults to 1
        :type verbose: int, optional
        """
        self.roi_size = roi_size
        self.verbose = verbose

        self.noverlap = int(noverlap)

        self.cent_dist_0 = self.roi_size[0] - self.noverlap
        self.cent_dist_1 = self.roi_size[1] - self.noverlap
        
        if video is not None:
            self.image = video.mraw[0]
            self.pick_window()
        else:
            print('set the polygon points in self.polygon and call the `get_roi_grid` method')

    def pick_window(self):
        # Tkinter root and matplotlib figure
        root = tk.Tk()
        root.title('Pick points')
        fig = Figure(figsize=(15, 7))
        ax = fig.add_subplot(111)
        ax.grid(False)
        ax.imshow(self.image, cmap='gray')
        plt.show()

        # Embed figure in tkinter winodw
        canvas = FigureCanvasTkAgg(fig, root)
        canvas.get_tk_widget().pack(side='top', fill='both', expand=1)
        NavigationToolbar2Tk(canvas, root)

        # Initiate polygon
        self.polygon = [[], []]
        line, = ax.plot(self.polygon[1], self.polygon[0], 'r.-')

        if self.verbose:
            print('SHIFT + LEFT mouse button to pick a pole.\nRIGHT mouse button to erase the last pick.')

        self.shift_is_held = False

        def on_key_press(event):
            """Function triggered on key press (shift)."""
            if event.key == 'shift':
                self.shift_is_held = True

        def on_key_release(event):
            """Function triggered on key release (shift)."""
            if event.key == 'shift':
                self.shift_is_held = False

        def onclick(event):
            if event.button == 1 and self.shift_is_held:
                if event.xdata is not None and event.ydata is not None:
                    self.polygon[1].append(int(np.round(event.xdata)))
                    self.polygon[0].append(int(np.round(event.ydata)))
                    if self.verbose:
                        print(f'y: {np.round(event.ydata):5.0f}, x: {np.round(event.xdata):5.0f}')

            elif event.button == 3 and self.shift_is_held:
                if self.verbose:
                    print('Deleted the last point...')
                del self.polygon[1][-1]
                del self.polygon[0][-1]

            line.set_xdata(self.polygon[1])
            line.set_ydata(self.polygon[0])
            fig.canvas.draw()

        def handle_close(event):
            """On closing."""
            self.polygon = np.asarray(self.polygon).T
            if self.verbose:
                for i, point in enumerate(self.polygon):
                    print(f'{i+1}. point: x ={point[1]:5.0f}, y ={point[0]:5.0f}')
            
            self.points = self.get_roi_grid()

        # Connecting functions to event manager
        fig.canvas.mpl_connect('key_press_event', on_key_press)
        fig.canvas.mpl_connect('key_release_event', on_key_release)
        fig.canvas.mpl_connect('button_press_event', onclick)
        # on closing the figure
        fig.canvas.mpl_connect('close_event', handle_close)

        root.mainloop()
    
    def get_roi_grid(self):
        points = self.polygon
        
        low_0 = np.min(points[:, 0])
        high_0 = np.max(points[:, 0])
        low_1 = np.min(points[:, 1])
        high_1 = np.max(points[:, 1])
        
        rois = []
        for i in range(low_0+self.cent_dist_0, high_0-self.cent_dist_0, self.cent_dist_0):
            for j in range(low_1+self.cent_dist_1, high_1-self.cent_dist_1, self.cent_dist_1):
                if inside_polygon(i, j, self.polygon):
                    rois.append([i, j])
        return np.asarray(rois)


class Select:
    def __init__(self, video=None, roi_size=(11, 11), noverlap=0):
        self.verbose = 0
        self.shift_is_held = False
        
        self.roi_size = roi_size
        self.noverlap = int(noverlap)
        self.cent_dist_0 = self.roi_size[0] - self.noverlap
        self.cent_dist_1 = self.roi_size[1] - self.noverlap

        self.polygon = [[], []]
        self.deselect_polygon = [[], []]

        root = tk.Tk()
        root.title('Selection')

        self.screen_width = root.winfo_screenwidth()
        self.screen_height = root.winfo_screenheight()
        root.geometry(f'{int(0.9*self.screen_width)}x{int(0.9*self.screen_height)}')

        self.options = SelectOptions(root, self)
        button1 = tk.Button(root, text='Options', command=lambda: self.open_options(root))
        button1.pack(side='top')

        self.fig = Figure(figsize=(10, 7))
        self.ax = self.fig.add_subplot(111)
        self.ax.grid(False)
        self.ax.imshow(video.mraw[0], cmap='gray')

        # Initiate polygon
        self.line, = self.ax.plot([], [], 'r.-')
        self.line2, = self.ax.plot([], [], 'bo')

        plt.show(block=False)

        # Embed figure in tkinter winodw
        canvas = FigureCanvasTkAgg(self.fig, root)
        canvas.get_tk_widget().pack(side='top', fill='both', expand=1, padx=5, pady=5)
        NavigationToolbar2Tk(canvas, root)
        
        if self.verbose:
            print('SHIFT + LEFT mouse button to pick a pole.\nRIGHT mouse button to erase the last pick.')

        # Connecting functions to event manager
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)
        self.fig.canvas.mpl_connect('key_release_event', self.on_key_release)

        self.update_variables()

        tk.mainloop()

    def _mode_selection_polygon(self, get_rois=True):
        """Select polygon to compute the points based on ROI size and
        ROI overlap."""
        def onclick(event):
            if event.button == 1 and self.shift_is_held:
                if event.xdata is not None and event.ydata is not None:
                    if self.polygon[0]:
                        del self.polygon[1][-1]
                        del self.polygon[0][-1]

                    self.polygon[1].append(int(np.round(event.xdata)))
                    self.polygon[0].append(int(np.round(event.ydata)))

                    if self.polygon[0]:
                        self.polygon[1].append(self.polygon[1][0])
                        self.polygon[0].append(self.polygon[0][0])
                        
                    if self.verbose:
                        print(f'y: {np.round(event.ydata):5.0f}, x: {np.round(event.xdata):5.0f}')

            elif event.button == 3 and self.shift_is_held:
                if self.verbose:
                    print('Deleted the last point...')
                del self.polygon[1][-1]
                del self.polygon[0][-1]

            self.line.set_xdata(self.polygon[1])
            self.line.set_ydata(self.polygon[0])
            self.fig.canvas.draw()

            if get_rois:
                self.plot_selection()
        
        self.fig.canvas.mpl_connect('button_press_event', onclick)

    def on_key_press(self, event):
        """Function triggered on key press (shift)."""
        if event.key == 'shift':
            self.shift_is_held = True
    
    def on_key_release(self, event):
        """Function triggered on key release (shift)."""
        if event.key == 'shift':
            self.shift_is_held = False

    def get_roi_grid(self, polygon_points):
        points = np.array(polygon_points)
        if points.shape[0] == 2:
            points = points.T
        
        low_0 = np.min(points[:, 0])
        high_0 = np.max(points[:, 0])
        low_1 = np.min(points[:, 1])
        high_1 = np.max(points[:, 1])
        
        rois = []
        for i in range(low_0+self.cent_dist_0, high_0, self.cent_dist_0):
            for j in range(low_1+self.cent_dist_1, high_1, self.cent_dist_1):
                if inside_polygon(i, j, points):
                    rois.append([i, j])
        return np.asarray(rois)
    
    def update_variables(self):
        self.line2.set_xdata([])
        self.line2.set_ydata([])
        self.fig.canvas.draw()

        self.mode = self.options.combobox.get()
        if self.mode == 'Selection polygon - ROI grid':
            self._mode_selection_polygon()

            self.roi_size = [int(self.options.roi_entry_y.get()), int(self.options.roi_entry_x.get())]
            self.noverlap = int(self.options.noverlap_entry.get())

            self.cent_dist_0 = self.roi_size[0] - self.noverlap
            self.cent_dist_1 = self.roi_size[1] - self.noverlap

            self.plot_selection()

        elif self.mode == 'Select polygon':
            self._mode_selection_polygon(get_rois=False)
        else:
            raise Exception('Non existing mode...')
        
    def plot_selection(self):
        if self.polygon[0] and self.polygon[1]:
            self.points = self.get_roi_grid(self.polygon)

            self.line2.set_xdata(self.points[:, 1])
            self.line2.set_ydata(self.points[:, 0])
            self.fig.canvas.draw()

    def clear_selection(self):
        self.polygon = [[], []]
        self.points = [[], []]
        self.clear_plot()
    
    def clear_plot(self):
        self.line.set_xdata([])
        self.line.set_ydata([])
        self.line2.set_xdata([])
        self.line2.set_ydata([])
        self.fig.canvas.draw()

    def open_options(self, root):
        self.options = SelectOptions(root, self)

class SelectOptions:
    def __init__(self, root, parent):
        root1 = tk.Toplevel(root)
        root1.title('Selection options')
        root1.geometry(f'{int(0.2*parent.screen_width)}x{int(0.5*parent.screen_height)}')

        roi_x = tk.StringVar(root1, value='21')
        roi_y = tk.StringVar(root1, value='21')
        noverlap = tk.StringVar(root1, value='0')

        row = 0
        tk.Label(root1, text='Selection mode:').grid(row=row, column=0)
        self.combobox = ttk.Combobox(root1, values = [
            'Select polygon', 
            'Selection polygon - ROI grid',
            ])
        self.combobox.current(0)
        self.combobox.grid(row=row, column=1, sticky='wens', padx=5, pady=5)

        row = 1
        tk.Label(root1, text='Horizontal ROI size').grid(row=row, column=0, sticky='E')
        self.roi_entry_x = tk.Entry(root1, textvariable=roi_x)
        self.roi_entry_x.grid(row=row, column=1, padx=5, pady=5, sticky='W')

        row = 2
        tk.Label(root1, text='Vertical ROI size').grid(row=row, column=0, sticky='E')
        self.roi_entry_y = tk.Entry(root1, textvariable=roi_y)
        self.roi_entry_y.grid(row=row, column=1, padx=5, pady=5, sticky='W')

        row = 3
        tk.Label(root1, text='Overlap pixels').grid(row=row, column=0, sticky='E')
        self.noverlap_entry = tk.Entry(root1, textvariable=noverlap)
        self.noverlap_entry.grid(row=row, column=1, padx=5, pady=5, sticky='W')

        # row = 4
        # self.canvas = tk.Canvas(root1, width=int(0.18*parent.screen_width))
        # self.canvas.grid(row=row, column=0, columnspan=2)

        # self.canvas.create_text(100, 100, text='Teset', font="Times 12")

        row = 4
        apply_button = tk.Button(root1, text='Apply', command=parent.update_variables)
        apply_button.grid(row=row, column=0, sticky='we', padx=5, pady=5)

        clear_button = tk.Button(root1, text='Clear', command=parent.clear_selection)
        clear_button.grid(row=row, column=1, sticky='w', padx=5, pady=5)


def inside_polygon(x, y, points):
    """
    Return True if a coordinate (x, y) is inside a polygon defined by
    a list of verticies [(x1, y1), (x2, x2), ... , (xN, yN)].

    Reference: http://www.ariel.com.au/a/python-point-int-poly.html
    """
    n = len(points)
    inside = False
    p1x, p1y = points[0]
    for i in range(1, n + 1):
        p2x, p2y = points[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / \
                            (p2y - p1y) + p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y
    return inside