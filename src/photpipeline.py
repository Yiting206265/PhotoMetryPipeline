from __future__ import print_function
import photutils
from photutils import DAOStarFinder
from astropy.stats import mad_std
from photutils import aperture_photometry, CircularAperture
from fitsimg import FitsImage
from astropy.stats import mad_std
import matplotlib.pyplot as plt
from astropy.visualization import LogStretch, SqrtStretch, AsinhStretch, HistEqStretch,ZScaleInterval
from astropy.visualization.mpl_normalize import ImageNormalize
import gkastro
from photutils import CircularAnnulus
import fitsimg
import numpy as np
import pandas as pd
import astropy
import os
from astroscrappy import detect_cosmics

class PhotoMetryImage(FitsImage):
    """
    A class to do photometry with fits images

    EXAMPLE:
    df = pd.read_csv("/Volumes/macwin/APORUN/APODATA_20161101/Q2OU02/UT170420/0_SCIENCE/TESTING_pipeline/RESULTS/aij_34_60_90_noPlane_noSigma_noRADEC2_measurements.txt",sep="\t")
    xx = df[["X(FITS)_T1","X(FITS)_C2","X(FITS)_C3","X(FITS)_C4","X(FITS)_C5","X(FITS)_C6","X(FITS)_C7"]].ix[0].values
    yy = df[["Y(FITS)_T1","Y(FITS)_C2","Y(FITS)_C3","Y(FITS)_C4","Y(FITS)_C5","Y(FITS)_C6","Y(FITS)_C7"]].ix[0].values
    PMI.perform_photometry_window_centroid(xx,yy,method="centroid_com")
    """
    saturation_warning_threshold = 55000.
    
    def __init__(self,filename):
        self.fimg = fitsimg.FitsImage(filename)
        self.data = self.fimg.data
        self.filename = filename
        self.basename = self.filename.split(os.sep)[-1]
        #self.bkg_sigma = mad_std(self.data)
                

    def remove_cosmics(self):
        self.cosmics_mask, self.data = detect_cosmics(self.data)
        self.num_cosmics = self.cosmics_mask.sum()
        print("Removing cosmics, found:", PMI.num_cosmics)

    #def calc_background(self):
    #    self.bkg_sigma = mad_std(self.data)
        
    def plot_apertures(self,stretch="hist",cmap="gray",alpha=0.5,lw=1.5,**kwargs):
        """
        Plot the apertures overlayed on the image
        """
        #self.fig, self.ax = plt.subplots()
        #self.ax.imshow(self.data, cmap='gray_r', origin='lower')
        self.plot(stretch=stretch,cmap=cmap)
        if type(self.apertures)==list:
            [ap.plot(color='blue',alpha=alpha,ax=self.ax,lw=1.5,**kwargs) for ap in self.apertures]
        else:
            self.apertures.plot(color='blue',alpha=alpha,ax=self.ax,lw=lw,**kwargs)

    def plot_photometry(self):
        """
        """
        self.fig, self.ax = plt.subplots()
        self.ax.plot(phot_table)

    def howell_center(self,postage_stamp):
        """
        Howell centroiding, from Howell's Handbook of CCD astronomy

        INPUT:
         postage_stamp - A 2d numpy array to do the centroiding

        OUTPUT:
         x and y center of the numpy array

        NOTES:
        Many thanks to Thomas Beatty and the MINERVAphot.py pipeline for this method
        see here: https://github.com/TGBeatty/MINERVAphot/blob/master/MINERVAphot.py
        """
        xpixels = np.arange(postage_stamp.shape[1])
        ypixels = np.arange(postage_stamp.shape[0])
        I = np.sum(postage_stamp, axis=0)
        J = np.sum(postage_stamp, axis=1)

        Isub = I-np.sum(I)/I.size
        Isub[Isub<0] = 0
        Jsub = J-np.sum(J)/J.size
        Jsub[Jsub<0] = 0
        xc = np.sum(Isub*xpixels)/np.sum(Isub)
        yc = np.sum(Jsub*ypixels)/np.sum(Jsub)
        return xc, yc
                
    def get_centroid_cutout(self,x,y,box_size=30,method="howell",dao_fwhm=10.,dao_SNR=100.,plot=False,plot_full=False,stretch=None):
        """
        Perform centroiding on a cutout window.
        
        INPUT:
            x -  a touple 
            y - 
            box_size -  the total box size in pixels
            method   -  the centroiding method to use for cutout. This can be
                            "centroid_2dg"
                            "centroid_com"
                            "daofind"
        OUTPUT:
            centroid positions
        
        EXAMPLE:
            PM = PhotoMetry(filenames[0])
            PM.perform_photometry()
            print(PM.get_centroid_cutout(234.227204448,25.6594515323,method="daofind",plot=True))
            print(PM.get_centroid_cutout(234.227204448,25.6594515323,method="centroid_com",plot=True))
            print(PM.get_centroid_cutout(234.227204448,25.6594515323,method="centroid_2dg",plot=True))

            print(PM.get_centroid_cutout(234.227204448,25.6594515323,method="daofind",box_size=50,stretch="HistEqStretch"))
            print(PM.get_centroid_cutout(234.227204448,25.6594515323,method="centroid_com",plot=True,box_size=50))
            print(PM.get_centroid_cutout(234.227204448,25.6594515323,method="centroid_2dg",plot=True,box_size=50))
        """
        int_x = int(round(x))
        int_y = int(round(y))
        postage_stamp = photutils.utils.cutouts.cutout_footprint(self.data,(int_x,int_y),box_size)
        
        # Only interested in the image data [1] is the mask
        postage_stamp = np.array(postage_stamp[0])
        #plt.imshow(postage_stamp)
        
        #print(postage_stamp)
        if method=="centroid_2dg":
            x_stamp_centroid, y_stamp_centroid = photutils.centroid_com(postage_stamp)
        elif method=="centroid_com":
            x_stamp_centroid, y_stamp_centroid = photutils.centroid_2dg(postage_stamp)
        elif method=="daofind":
            daofind = DAOStarFinder(fwhm = dao_fwhm, threshold=dao_SNR*self.bkg_sigma)
            sources = daofind(postage_stamp)
            positions = (sources['xcentroid'], sources['ycentroid']) 
            x_stamp_centroid, y_stamp_centroid = float(positions[0]),float(positions[1])
        elif method=="howell":
            x_stamp_centroid, y_stamp_centroid = self.howell_center(postage_stamp)
        else:
            print("Error: method must be 'daofind', centroid_2dg', 'centroid_com' or 'howell'")
            pass
        x_centroid = x_stamp_centroid + int_x - box_size/2.
        y_centroid = y_stamp_centroid + int_y - box_size/2.
        if plot:
            fig, ax = plt.subplots()
            if stretch:
                norm = gkastro.stretch_data(postage_stamp,method=stretch)
            else:
                norm = None
            ax.imshow(postage_stamp,origin="lower",extent=[int_x-box_size/2.,
                                                           int_x+box_size/2.,
                                                           int_y-box_size/2.,
                                                           int_y+box_size/2.],interpolation="none",norm=norm)
            ax.set_xlim(int_x-box_size/2.,int_x+box_size/2.)
            ax.set_ylim(int_y-box_size/2.,int_y+box_size/2.)
            ax.set_xlabel("X pixel")
            ax.set_ylabel("Y pixel")
            ax.plot(x_centroid, y_centroid, color='#1f77b4', marker="+", ms=30, mew=2)
        if plot_full:
            fig, ax = plt.subplots()
            ax.imshow(self.data,origin="lower")
            ax.plot(x_centroid, y_centroid, color='#1f77b4', marker="+", ms=30, mew=2)
        return x_centroid, y_centroid
        
    def perform_photometry_window_centroid(self,xcen,ycen,r_aper=34.,r_annulus1=60.,r_annulus2=90.,
                                           box_size=80.,method="howell"):
        """
        Perform the photometry by window centroiding

        Currently only supports one aperture radius
        """
        if (r_annulus1 < r_aper) ^ (r_annulus2 < r_aper):
            raise Exception, '(r_annulus1 < r_aper) ^ (r_annulus2 < r_aper)'
            
        box_size = int(box_size)
        self.new_xpos = np.zeros(len(xcen))
        self.new_ypos = np.zeros(len(ycen))
        
        # Recentroid on each of the reference star
        for i in range(len(xcen)):
            self.new_xpos[i], self.new_ypos[i] = self.get_centroid_cutout(xcen[i],ycen[i],box_size=box_size,method=method,plot=False,plot_full=False)
        
        # Save new centroids
        self.positions = np.vstack([self.new_xpos,self.new_ypos]).T
        
        # Create apertures
        self.apertures = CircularAperture2(self.positions,r=r_aper)
        self.annulus_apertures = CircularAnnulus2(self.positions, r_in=r_annulus1,r_out=r_annulus2)
        apers = [self.apertures, self.annulus_apertures]
        
        # Perform photometry
        self.phot_table = aperture_photometry(self.data,apers,method="exact",subpixels=5)
        
        self.phot_table["Peak"]   = self.apertures.get_peak_in_aperture(self.data)
        self.data_cutouts_aper    = self.apertures.get_data_cutouts(self.data)
        self.data_cutouts_annulus = self.annulus_apertures.get_data_cutouts(self.data)
        
        # 
        bkg_mean  = self.phot_table['aperture_sum_1'] / self.annulus_apertures.area()
        bkg_sum   = bkg_mean * self.apertures.area()
        final_sum = self.phot_table['aperture_sum_0'] - bkg_sum
        
        self.phot_table["final_sum"] = final_sum
        self.phot_table["bkg_mean"]  = bkg_mean
        
        # master dataframe
        df = pd.DataFrame(index=[1]) # only one row

        df["Label"]                    = self.fimg.header["FILENAME"]
        df["JD_UTC"]                   = astropy.time.Time(self.fimg.header["DATE-OBS"],format="isot").jd + (self.fimg.header["EXPTIME"]/2.)/(24.*3600.)
        df["JD_SOBS"]                  = astropy.time.Time(self.fimg.header["DATE-OBS"],format="isot").jd
        df["AIRMASS"]                  = self.fimg.header["AIRMASS"]
        df["TELAZ"]                    = self.fimg.header["TELAZ"]
        df["TELALT"]                   = self.fimg.header["TELAZ"]
        df["TELROT"]                   = self.fimg.header["TELROT"]
        
        # Loop over the stars
        for i in range(len(xcen)):
            if i==0:
                TorC = "T"
            else:
                TorC = "C"
            df["X(FITS)_"+TorC+str(i+1)]   = self.phot_table["xcenter"][i].value
            df["Y(FITS)_"+TorC+str(i+1)]   = self.phot_table["ycenter"][i].value
            df["Sky/Pixel_"+TorC+str(i+1)] = self.phot_table["bkg_mean"][i]
            df["Source-Sky_"+TorC+str(i+1)]= self.phot_table["final_sum"][i]
            df["Peak_"+TorC+str(i+1)]      = self.phot_table["Peak"][i]
        
        self.source_sky_C_keys = [key for key in df.keys() if "Source-Sky_C" in key]
        self.sky_pixel_keys    = [key for key in df.keys() if "Sky/Pixel" in key]
        df["tot_C_cnts"] = df[self.source_sky_C_keys].sum().sum()
        df["rel_flux_T1"]= df["Source-Sky_T1"]/df["tot_C_cnts"]
        
        for i in range(len(xcen)):
            if i==0:
                TorC = "T"
            else:
                TorC = "C"
            df["rel_flux_"+TorC+str(i+1)]=df["Source-Sky_"+TorC+str(i+1)]/df["tot_C_cnts"]
            
        peak_flux = np.max(df[[key for key in df.keys() if "Peak" in key]].values)
        if peak_flux > self.saturation_warning_threshold:
            df["Saturated"] = peak_flux
        else:
            df["Saturated"] = 0.
        return df
    
    
class PhotoMetryPipeLine(object):
    """
    A class to perform photometry, looping through PhotoMetryImage classes

    # STEP 1:
      - Find reference stars using AIJ. Save file

    # STEP 2:
    df = pd.read_csv("/Volumes/macwin/APORUN/APODATA_20161101/Q2OU02/UT170420/0_SCIENCE/TESTING_pipeline/RESULTS/aij_34_60_90_noPlane_noSigma_noRADEC2_measurements.txt",sep="\t")
    xx = df[["X(FITS)_T1","X(FITS)_C2","X(FITS)_C3","X(FITS)_C4","X(FITS)_C5","X(FITS)_C6","X(FITS)_C7"]].ix[0].values
    yy = df[["Y(FITS)_T1","Y(FITS)_C2","Y(FITS)_C3","Y(FITS)_C4","Y(FITS)_C5","Y(FITS)_C6","Y(FITS)_C7"]].ix[0].values
    PMI.perform_photometry_window_centroid(xx,yy,method="centroid_com")
    """
    def __init__(self,filenames,xcenters,ycenters,r_aper=34.,r_annulus1=60.,r_annulus2=90.,box_size=80.,method="centroid_com"):
        #self.foldername = foldername
        #self.regex = regex
        self.filenames = filenames
        self.numfiles = len(self.filenames)
        self.xcenters = xcenters
        self.ycenters = ycenters
        self.r_aper = r_aper
        self.r_annulus1 = r_annulus1
        self.r_annulus2 = r_annulus2
        self.method = method
        self.box_size = box_size
        
    def run_pipeline(self,remove_cosmics=False):
        """
        Run the pipeline

        TODO:
        - Fix initilization of centroids
        - Add plotting capabilities
        """
        df_index = range(self.numfiles)

        # Loop over files
        for i,filename in enumerate(self.filenames):
            PMI = PhotoMetryImage(self.filenames[i])
            print("Analyzing file #",i,PMI.basename,"with r=",self.r_aper,self.r_annulus1,self.r_annulus2)

            # This can only be done before calibration
            if remove_cosmics:
                PMI.remove_cosmics()
                        
            if i==0:
                self.df_master = PMI.perform_photometry_window_centroid(self.xcenters,self.ycenters,r_aper=self.r_aper,r_annulus1=self.r_annulus1,r_annulus2=self.r_annulus2,box_size=self.box_size,method=self.method)
                self.df_master = self.df_master.reset_index(drop=True)
            else:
                ## Perform windowed centroiding, using previous centroid
                #xkeys = [key for key in self.df_master if "X(FITS)_" in key]
                #ykeys = [key for key in self.df_master if "Y(FITS)_" in key]
                #xcenter = self.df_master[xkeys].ix[i-1].values  
                #ycenter = self.df_master[ykeys].ix[i-1].values  
                df_temp = PMI.perform_photometry_window_centroid(self.xcenters,self.ycenters,r_aper=self.r_aper,r_annulus1=self.r_annulus1,r_annulus2=self.r_annulus2,box_size=self.box_size,method=self.method)
                #df_temp = PMI.perform_photometry_window_centroid(xcenter,ycenter,r_annulus1=self.r_annulus1,r_annulus2=self.r_annulus2,box_size=self.box_size,method=self.method)
                self.df_master = pd.concat([self.df_master,df_temp])
                self.df_master = self.df_master.reset_index(drop=True)
        
        return self.df_master

    def save_photometry(self,outputname=""):
        if outputname=="":
            outputname = "phot_results_"+str(self.r_aper)+"_"+str(self.r_annulus1)+"_"+str(self.r_annulus2)+".csv"
        self.df_master.to_csv(outputname,index=False)
        print("Saved file to:",outputname)
    
      
    def rel_flux_residuals(self):
        df_pos = pd.read_csv("../data/Processed_science_2x2_binning_sdssi_fast_diff_measurements.xls",sep="\t")
        df_pos=df_pos[0:200]
        time=self.df_master[["JD_UTC"]]
        self.fig, self.ax = plt.subplots(nrows=3,sharex=True)#,figsize=(8,8))
        list_name=['rel_flux_T1_AIJ','rel_flux_T1_Pipline','Residuals']
        for i in np.arange(0,3):
            if i==0:
                normalized_flux_AIJ=df_pos["rel_flux_T1"]/np.nanmean(df_pos["rel_flux_T1"])
                self.ax[i].scatter(time,normalized_flux_AIJ,marker='.',color='r',label='AIJ result')
                self.ax[i].legend(prop={'size':8})
                self.ax[i].set_title("Relative flux vs. Time")
            if i==1:
                normalized_flux_Pipeline=self.df_master['rel_flux_T1']/np.nanmean(self.df_master['rel_flux_T1'])
                self.ax[i].scatter(time,normalized_flux_Pipeline,marker='.',color='b',label='Pipeline result')
                self.ax[i].legend(prop={'size':8})
            if i==2:
                residuals=df_pos["rel_flux_T1"]/np.nanmean(df_pos["rel_flux_T1"])-self.df_master['rel_flux_T1']/np.nanmean(self.df_master['rel_flux_T1'])
                self.ax[i].scatter(time,residuals,marker='.',color='g',label='Residuals: AIJ - Pipeline')
                self.ax[i].legend(prop={'size':8})
            self.ax[i].set_ylabel(list_name[i])
        #plt.savefig("./",format='png')    
         
            
    def plot_X_and_Y(self):
        #should change this to a double for-loop
        df_pos = pd.read_csv("../data/Processed_science_2x2_binning_sdssi_fast_diff_measurements.xls",sep="\t")
        df_pos=df_pos[0:200]
        time=self.df_master[["JD_UTC"]]
        self.fig, self.ax = plt.subplots(nrows=2,sharex=True)
        list_cen=["X(FITS)_T1","Y(FITS)_T1"]
        """
        for i in np.arange(0,2):
            if i==0:
                x_cen=self.df_master[list_cen[i]]
                self.ax[i].scatter(time,x_cen,marker='.',color='g',label='X centroid')
                self.ax[i].legend(prop={'size':8})
                self.ax[i].set_title("X and Y Centroids from Pipeline")
            if i==1:
                y_cen=self.df_master[list_cen[i]]
                self.ax[i].scatter(time,y_cen,marker='.',color='b',label='Y centroid')
                self.ax[i].legend(prop={'size':8})
            self.ax[i].set_ylabel(list_cen[i])
        plt.savefig("./centroids_pipeline.png",format='png')
        
        for i in np.arange(0,2):
            if i==0:
                x_cen=df_pos[list_cen[i]]
                self.ax[i].scatter(time,x_cen,marker='.',color='g',label='X centroid')
                self.ax[i].legend(prop={'size':8})
                self.ax[i].set_title("X and Y Centroids from AIJ")
            if i==1:
                y_cen=df_pos[list_cen[i]]
                self.ax[i].scatter(time,y_cen,marker='.',color='b',label='Y centroid')
                self.ax[i].legend(prop={'size':8})
            self.ax[i].set_ylabel(list_cen[i])
        plt.savefig("./centroids_AIJ.png",format='png')
        """
        for i in np.arange(0,2):
            if i==0:
                x_cen=df_pos[list_cen[i]]-self.df_master[list_cen[i]]
                self.ax[i].scatter(time,x_cen,marker='.',color='g',label='X centroid')
                self.ax[i].legend(prop={'size':8})
                self.ax[i].set_title("X and Y Residuals")
            if i==1:
                y_cen=df_pos[list_cen[i]]-self.df_master[list_cen[i]]
                self.ax[i].scatter(time,y_cen,marker='.',color='b',label='Y centroid')
                self.ax[i].legend(prop={'size':8})
            self.ax[i].set_ylabel(list_cen[i])
        plt.savefig("./centroids_difference.png",format='png')
        
    def plot_photometry(self):
        time=self.df_master[["JD_UTC"]]
        list_rel_flux = ["rel_flux_T1","rel_flux_C2","rel_flux_C3","rel_flux_C4","rel_flux_C5","rel_flux_C6"]
        #plot the raw photometry of the target star and reference stars
        num_refstars=5
        nrows=num_refstars+1
        self.fig, self.ax = plt.subplots(nrows=nrows,sharex=True,figsize=(8,8))
        
        for i in np.arange(0,nrows):
            if i==0:
                self.ax[i].scatter(time,self.df_master[list_rel_flux[i]]/np.nanmean(self.df_master[list_rel_flux[i]]),marker='.',color='r')
                self.ax[i].set_title("Photometry")
            else:
                self.ax[i].scatter(time,self.df_master[list_rel_flux[i]]/np.nanmean(self.df_master[list_rel_flux[i]]),marker='.')
            self.ax[i].set_ylabel(list_rel_flux[i])
        self.ax[num_refstars].set_xlabel("JD_UTC")
        #self.ax[i].ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        self.fig.tight_layout()


        
        


################## 
################## 
################## 
# Extending Circular Apertures
################## 
################## 

class CircularAperture2(CircularAperture):
    """
    Modifying so I can get the peak counts in the aperture
    """ 
    def __init__(self,positions,r):
        CircularAperture.__init__(self,positions,r)
        
    def get_data_cutouts(self,data,mask=None,method='exact',subpixels=5,unit=None):
        data = np.asanyarray(data)
        if mask is not None:
            mask = np.asanyarray(mask)
            data = copy.deepcopy(data)    # do not modify input data
            data[mask] = 0
        data_cutouts = []
        for mask in self.to_mask(method=method, subpixels=subpixels):
            data_cutout = mask.cutout(data)
            if data_cutout is None:
                data_cutouts.append(np.nan)
            else:
                data_cutouts.append(data_cutout * mask.data)
        return data_cutouts
    
    def get_peak_in_aperture(self,data,mask=None, method='exact',subpixels=5,unit=None):
        """
        Function to get the peak value in an aperture
        
        TODO:
        - Do this also for annulus ?
        """
        data_cutouts = self.get_data_cutouts(data,mask=None,method='exact',subpixels=5,unit=None)
        
        aperture_peaks = []
        for data_cutout in data_cutouts:
            aperture_peaks.append(np.max(data_cutout))
        
        aperture_peaks = self._prepare_photometry_output(aperture_peaks,
                                                        unit=unit)
        return aperture_peaks
    
class CircularAnnulus2(CircularAnnulus):
    
    def __init__(self,positions,r_in,r_out):
        CircularAnnulus.__init__(self,positions,r_in,r_out)
        
    def get_data_cutouts(self,data,mask=None,method='exact',subpixels=5,unit=None):
        data = np.asanyarray(data)
        if mask is not None:
            mask = np.asanyarray(mask)
            data = copy.deepcopy(data)    # do not modify input data
            data[mask] = 0
        data_cutouts = []
        for mask in self.to_mask(method=method, subpixels=subpixels):
            data_cutout = mask.cutout(data)
            if data_cutout is None:
                data_cutouts.append(np.nan)
            else:
                data_cutouts.append(data_cutout * mask.data)
        return data_cutouts
    
    def get_mean_background(self,sigma=3.):
        data_cutouts = self.get_data_cutouts(data,mask=None,method='exact',subpixels=5,unit=None)
        
        mean_backgrounds = []
        for data_cutout in data_cutouts:
            data = astropy.stats.sigma_clip(data_cutout,sigma=sigma)
            mean_backgrounds.append(data)
            
        
