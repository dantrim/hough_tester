#!/usr/vin/env python

import ROOT as r
r.gROOT.SetBatch(True)
r.gStyle.SetOptStat(False)

import numpy as np
import matplotlib.pyplot as plt

from math import cos, sin, sqrt

import sys

class Cluster :
    def __init__(self, x_, y_) :
        self._x = x_
        self._y = y_

    def x(self) :
        return self._x

    def y(self) :
        return self._y

    def Print(self) :
        print "Cluster   (%.2f, %.2f)"%(self._x, self._y)

def clusters_in_chamber(xc, yc) :
    x_out = []
    y_out = []

    for i, y in enumerate(yc) :
        if y>0 and y<5 :
            x_out.append(xc[i])
            y_out.append(y)

    return x_out, y_out

def make_clusters() :
    # offset_0
    offset_0 = 38
    slope_0 = 0.12

    np.random.seed(5)
    x = np.linspace(0,100,30)
    y0 = slope_0 * (x - offset_0) + np.random.normal(0,0.005,30)
    x_tmp, y_tmp = clusters_in_chamber(x, y0)

    out_clusters = []
    for i in xrange(len(x_tmp)) :
        c = Cluster(x_tmp[i], y_tmp[i])
        out_clusters.append(c)

    print "make_clusters    %d clusters"%len(out_clusters)
    return out_clusters



def draw_chamber_tracks(clusters) :

    c = r.TCanvas("c_cl","",1600,600)
    c.cd()
    hax = r.TH1F("hax","", 100,0,100)#,30,0,5)
    hax.GetYaxis().SetRangeUser(0,5)
    hax.Draw("axis")
    g = r.TGraph(0)
    for cl in clusters :
        x = cl.x()
        y = cl.y()
        print "%d  %.2f   %.2f"%(g.GetN(), x, y)
        g.SetPoint(g.GetN(), float(x), float(y))#, 1)
    g.SetMarkerStyle(20)
    g.SetMarkerColor(r.kBlue)
    g.Draw("p0")
    c.SaveAs("chamber_tracks.eps")

def draw_accumulator(accumulator) :
    
    fig = plt.figure()
    plt.imshow(accumulator, origin="lower", cmap="jet", interpolation="nearest", aspect="auto")
    fig.savefig("accumulator.eps")

def build_accumulator(clusters) :

    x_values = [cl.x() for cl in clusters]
    y_values = [cl.y() for cl in clusters]
    points = np.zeros((101,101))
    n_max = max(len(x_values), len(y_values))
    for i in xrange(n_max) :
        x_i = x_values[i]
        y_i = y_values[i]
        points[x_i, y_i] += 1

    thetas = np.deg2rad(np.arange(-90.0,90.0))
    width, height = points.shape
    diag_len = np.ceil(np.sqrt(width*width + height*height))
    #print "diag ",diag_len
    #rhos = np.linspace(-diag_len, diag_len, 2.0*diag_len)
    rhos = np.linspace(0,2*diag_len,2.0*diag_len)

    # cache the cosines and sines over theta range
    cos_t = np.cos(thetas)
    sin_t = np.sin(thetas)
    n_thetas = len(thetas)

    accumulator = np.zeros((2*diag_len, n_thetas), dtype=np.uint64)
    y_idxs, x_idxs = np.nonzero(points)

    for i in range(len(x_idxs)) :
        x = x_idxs[i]
        y = y_idxs[i]
        for t_idx in range(n_thetas) :
            rho = round(x*cos_t[t_idx] + y*sin_t[t_idx]) + diag_len
            accumulator[rho, t_idx] += 1

    draw_accumulator(accumulator)

    return accumulator, thetas, rhos

def get_cluster_group(theta_max, rho_max, thetas, rhos, clusters) :

    print "get_cluster_group"

    th_window = 5 #degrees
    rh_window = 10

    
    x_values = [cl.x() for cl in clusters]
    y_values = [cl.y() for cl in clusters]
    points = np.zeros((101,101))
    n_max = max(len(x_values), len(y_values))
    for i in xrange(n_max) :
        x_i = x_values[i]
        y_i = y_values[i]
        points[x_i, y_i] += 1

    width, height = points.shape
    diag_len = np.ceil(np.sqrt(width*width + height*height))

    # cache the cosines and sines over theta range
    cos_t = np.cos(thetas)
    sin_t = np.sin(thetas)
    n_thetas = len(thetas)


    ok_thetas = []
    for th in thetas :
        th = np.rad2deg(th)
        if th > (np.rad2deg(theta_max)-th_window) and th < (np.rad2deg(theta_max)+th_window) :
            ok_thetas.append(th)
    ok_rhos = []
    for rh in rhos :
        if rh > (rho_max-rh_window) and rh < (rho_max+rh_window) :
            ok_rhos.append(rh)

    print ok_thetas


    out_clusters = []
    idx_to_remove = []
    for icluster, cl in enumerate(clusters) :
        # get the hough line for this cluster
        xcl = cl.x()
        ycl = cl.y()
        for th in ok_thetas :
            th = np.deg2rad(th)
            print th
            rho_check = round(xcl*cos_t[th] + ycl*sin_t[th]) + diag_len
            print rho_check
            if rho_check > (rho_max-rh_window) and rh < (rho_max+rh_window) :
                print "INSIDE WINDOW"
                out_clusters.append(cl)
                idx_to_remove.append(icluster)
    print "rho window: %.2f < %.2f < %.2f"%(rho_max-rh_window, rho_max, rho_max+rh_window)
    print "the window: %.2f < %.2f < %.2f"%(np.rad2deg(theta_max)-th_window, np.rad2deg(theta_max), np.rad2deg(theta_max)+th_window)
    sys.exit()
    tmp_cl = []
    for icl, cl in enumerate(clusters) :
        if icl in idx_to_remove : continue
        tmp_cl.append(cl)
    clusters = tmp_cl

    return out_clusters, clusters
        
        
    
    

def main() :

    # get the clusters
    clusters = make_clusters()
    for cl in clusters :
        cl.Print()
    
    # draw the chamber tracks
    draw_chamber_tracks(clusters)

    # fill the accumulator
    accumulator, thetas, rhos = build_accumulator(clusters)

    # get the track-associated clusters
    track_cluster_groups = []
    while True :
        thetas = [th+(3.14/2) for th in thetas] 
        print thetas
        print rhos

        # get the area of maximum occupancy
        max_idx = np.argmax(accumulator)
        rho = rhos[max_idx / accumulator.shape[1]]
        theta = thetas[max_idx % accumulator.shape[1]]

        print "%d : (theta, rho) = (%.2f,%.2f)"%(max_idx, np.rad2deg(theta), rho)

        track_clusters, clusters = get_cluster_group(theta, rho, thetas, rhos, clusters)
        print len(track_clusters)
        track_cluster_groups.append(track_clusters)
        sys.exit()


#####################################################################
if __name__=="__main__" :
    main()
