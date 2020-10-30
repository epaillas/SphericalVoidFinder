import numpy as np
import sys
import os
import glob
import subprocess
from astropy.io import fits
from scipy.spatial import Delaunay
from scipy.io import FortranFile
import matplotlib.pyplot as plt


class SphericalVoids:

    def __init__(self,
                 centres_filename,
                 tracers_filename,
                 output_handle,
                 ncores,
                 steps,
                 pos_cols,
                 rvoidmax,
                 box_size,
                 ngrid):

        steps = [int(i) for i in steps.split(',')]
        pos_cols = [int(i) for i in pos_cols.split(',')]

        self.centres_filename = centres_filename
        self.tracers_filename = tracers_filename
        self.output_handle = output_handle
        self.ncores = ncores
        self.steps = steps
        self.pos_cols = pos_cols
        self.rvoidmax = rvoidmax
        self.box_size = box_size
        self.density_threshold = 0.2
        self.ngrid = ngrid

        if 1 in steps:
            self.get_circumcentres()

        # Grow spheres from the centres found in the previous step
        if 2 in steps:
            self.grow_spheres()

        # Find better void centres by shifting the original positions
        if 3 in steps:
            self.recentre_spheres()

        # Sort spheres in decreasing order of radius
        # and remove overlapping spheres
        if 4 in steps:
            self.sort_spheres()
            self.overlap_filter(overlap=0.0)
            self.overlap_filter(overlap=0.2)
            self.overlap_filter(overlap=0.5)

    def concat_files(self, input_files, output_file):
        with open(output_file, 'w+') as outfile:
            for fname in input_files:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

    def get_periodic_images(self, data):
        '''
        Find the relevant images of a 
        set of points in a box that
        has boundary conditions.
        '''
        images = []
        buffer = self.box_size / 10  # probably an overkill

        for point in data:
            condx = ((self.box_size - buffer) <
                     point[0]) or (point[0] < buffer)
            condy = ((self.box_size - buffer) <
                     point[1]) or (point[1] < buffer)
            condz = ((self.box_size - buffer) <
                     point[2]) or (point[2] < buffer)

            if condx and condy and condz:
                shiftx = point[0] + \
                    np.copysign(self.box_size, (buffer - point[0]))
                shifty = point[1] + \
                    np.copysign(self.box_size, (buffer - point[1]))
                shiftz = point[2] + \
                    np.copysign(self.box_size, (buffer - point[2]))
                images.append([shiftx, shifty, shiftz])
            if condx and condy:
                shiftx = point[0] + \
                    np.copysign(self.box_size, (buffer - point[0]))
                shifty = point[1] + \
                    np.copysign(self.box_size, (buffer - point[1]))
                shiftz = point[2]
                images.append([shiftx, shifty, shiftz])
            if condx and condz:
                shiftx = point[0] + \
                    np.copysign(self.box_size, (buffer - point[0]))
                shifty = point[1]
                shiftz = point[2] + \
                    np.copysign(self.box_size, (buffer - point[2]))
                images.append([shiftx, shifty, shiftz])
            if condy and condz:
                shiftx = point[0]
                shifty = point[1] + \
                    np.copysign(self.box_size, (buffer - point[1]))
                shiftz = point[2] + \
                    np.copysign(self.box_size, (buffer - point[2]))
                images.append([shiftx, shifty, shiftz])
            if condx:
                shiftx = point[0] + \
                    np.copysign(self.box_size, (buffer - point[0]))
                shifty = point[1]
                shiftz = point[2]
                images.append([shiftx, shifty, shiftz])
            if condy:
                shiftx = point[0]
                shifty = point[1] + \
                    np.copysign(self.box_size, (buffer - point[1]))
                shiftz = point[2]
                images.append([shiftx, shifty, shiftz])
            if condz:
                shiftx = point[0]
                shifty = point[1]
                shiftz = point[2] + \
                    np.copysign(self.box_size, (buffer - point[2]))
                images.append([shiftx, shifty, shiftz])

        images = np.asarray(images)
        return images

    def get_circumcentres(self, radius_limit=1000):
        '''
        Find the centre of the circumspheres
        associated to an input catalogue of
        tetrahedra.
        '''

        print('Finding circumcentres of tetrahedra...')
        vertices = self.delaunay_triangulation()
        cenx, ceny, cenz, r = [], [], [], []

        sing = 0
        for tetra in vertices:
            x0, x1, x2, x3 = tetra
            A = []
            B = []
            A.append((x1 - x0).T)
            A.append((x2 - x0).T)
            A.append((x3 - x0).T)
            A = np.asarray(A)
            B = np.sum(A**2, axis=1)
            B = np.asarray(B)
            try:
                C = np.linalg.inv(A).dot(B)
            except:
                sing += 1
                continue
            centre = x0 + 0.5 * C
            radius = 0.5 * np.sqrt(np.sum(C**2))
            if radius < radius_limit:
                cenx.append(centre[0])
                ceny.append(centre[1])
                cenz.append(centre[2])
                r.append(radius)

        print('{} singular matrices found.'.format(sing))

        cenx = np.asarray(cenx)
        ceny = np.asarray(ceny)
        cenz = np.asarray(cenz)

        cenx = cenx.reshape(len(cenx), 1)
        ceny = ceny.reshape(len(ceny), 1)
        cenz = cenz.reshape(len(cenz), 1)

        cout = np.hstack([cenx, ceny, cenz])

        print('{} centres found.'.format(len(cout)))

        f = FortranFile(self.centres_file, 'w')
        nrows, ncols = np.shape(cout)
        f.write_record(nrows)
        f.write_record(ncols)
        f.write_record(cout)
        f.close()

        self.centres = cout

        return

    def delaunay_triangulation(self, guards=False):
        '''
        Make a Delaunay triangulation over
        the cartesian positions of the tracers.
        Returns the vertices of tetrahedra.
        '''
        f = FortranFile(self.tracers_filename, 'r')
        nrows = f.read_ints()[0]
        ncols = f.read_ints()[0]
        data = f.read_reals(dtype=np.float64).reshape(nrows, ncols)
        f.close()

        points = data[:, :3]

        # add periodic images
        images = self.get_periodic_images(points)
        points = np.vstack([points, images])

        triangulation = Delaunay(points)
        simplices = triangulation.simplices.copy()
        vertices = points[simplices]
        print('{} vertices found.'.format(len(vertices)))
        return vertices

    def grow_spheres(self):
        '''
        Grow spheres from an input 
        catalogue of centres.
        '''
        print('Proceeding to grow spheres...')
        binpath = sys.path[0] + '/bin/'
        output_filename = self.output_handle + '.SVF'

        cmd = ['mpirun',
               '-np',
               str(self.ncores),
               binpath + 'grow_spheres.exe',
               self.tracers_filename,
               self.centres_filename,
               output_filename,
               str(self.box_size),
               str(self.density_threshold),
               str(self.rvoidmax),
               str(self.ngrid)]

        logfile = self.output_handle + '_grow_spheres.log'
        log = open(logfile, "w+")

        subprocess.call(cmd, stdout=log, stderr=log)

        if self.ncores > 1:
            files = glob.glob(output_filename + '.*')
            self.concat_files(input_files=files, output_file=output_filename)
            subprocess.call(['rm'] + files)

        return

    def recentre_spheres(self):
        '''
        Find better centres for an input
        void catalogue.
        '''
        print('Recentring spheres...')
        binpath = sys.path[0] + '/bin/'
        centres_filename = self.output_handle + '.SVF'
        output_filename = self.output_handle + '.SVF_recen'

        cmd = ['mpirun',
               '-np',
               str(self.ncores),
               binpath + 'recentring.exe',
               self.tracers_filename,
               centres_filename,
               output_filename,
               str(self.box_size),
               str(self.density_threshold),
               str(self.rvoidmax),
               str(self.ngrid)]

        logfile = self.output_handle + '_recentring.log'
        log = open(logfile, "w+")

        subprocess.call(cmd, stdout=log, stderr=log)

        if self.ncores > 1:
            files = glob.glob(output_filename + '.*')
            self.concat_files(input_files=files,
                              output_file=output_filename)
            subprocess.call(['rm'] + files)

        return

    def sort_spheres(self, fname='', radius_col=3):
        '''
        Sort an input void catalogue in
        decreasing order of radius.
        '''
        print('Sorting spheres by decreasing radius...')
        if fname == '':
            fname = self.output_handle + '.SVF_recen'

        voids = np.genfromtxt(fname)
        voids = voids[np.argsort(voids[:, radius_col])]
        voids = voids[::-1]

        fmt = 4*'%10.3f ' + '%10i ' + '%10.3f '
        np.savetxt(fname, voids, fmt=fmt)

        return

    def overlap_filter(self, overlap=0.0):

        binpath = sys.path[0] + '/bin/'

        input_filename = self.output_handle + '.SVF_recen'
        output_filename = self.output_handle + \
            'SVF_recen_ovl{}'.format(overlap)

        cmd = [binpath + 'overlapping.exe',
               input_filename,
               output_filename,
               str(self.box_size),
               str(overlap),
               str(self.ngrid)]

        logfile = self.output_handle + '_overlapping.log'
        log = open(logfile, "w+")

        subprocess.call(cmd, stdout=log, stderr=log)

        return
