#!/usr/bin/env python
from tools import rdarg
from sys import argv
from numpy import arange
try:
    from lib.pp import pp as ppserver
except:
    pass
import time
from copy import copy,deepcopy
import gc
import thread
from science import science
import os
import socket
try:
    import matplotlib
    import matplotlib.pyplot as plt
except:
    print 'Could not init. graphics - did you log in with ssh -X? (or Y)?'
import pickle
from dataIO import getData,pp_run,save,load,Progress,sel_edit,sel_edit2,pp_kill
from dataIO import getcfgProp, setcfgProp

import os
import numpy
from math import pi
from detections import accelerate

class getOut(Exception):
    import os

    def __init__(self,msg):
        print '\t#########################################################################'
        print '\t\n\n\n\n\t\t%s\n\n\n\n'%(msg)
        print '\t#########################################################################'
#        raise SystemExit(0)




####################################################################################
############################## PP DEPENDENT FUNCTIONS ##############################
####################################################################################
from parallel_detect import parallel_detect
from parallel_vertices import deploy_vertices
from parallel_qhull import qhull
from science import cherrypick
import dataIO

if 1:
    global range,intercept,colours,gradient,ncpu,dc,width,co_dets,init_width,fileIO,sci,cores,all_collapsed_clusters,refinery,job_server,source,export,zoo,profile,msgs,subset_index,to_disk,loc,version_id,codets_disk,cluster_zoo,gladders_fit,metric,range,stripe,do_train,cosma,magmod,dmag,cg_thresh,pe_thresh,train_path,do_evo,is_mock

def run(data,core,iteration):
    global range,intercept,colours,gradient,ncpu,dc,width,co_dets,init_width,fileIO,sci,cores,all_collapsed_clusters,refinery,job_server,source,export,zoo,profile,msgs,subset_index,to_disk,loc,version_id,codets_disk,cluster_zoo,gladders_fit,metric,range,stripe,do_train,cosma,magmod,dmag,cg_thresh,pe_thresh,train_path,do_evo,is_mock
    import numpy
    ####################################################################################
    ############################## PARALLEL PYTHON INIT ################################
    ####################################################################################
    ncpus = ncpu
    pp = True
    
    do_refine = False
    if ncpus == 0 or ncpus ==1: pp = False
    print 'Parallel-python: %s'%(str(pp))
    save_first_stage = True


    completed = []
    #ppservers=("*",)              # use for autodiscovery of clusters nodes
    ppservers = ()
    dependent_functions = (deploy_vertices,qhull,parallel_detect,sci.richness,cherrypick,dataIO.Progress,sel_edit2)
    modules = ('science','parallel_detect','clusters','numpy','parallel_qhull','parallel_vertices','vperc','time','detections')
    source = {}
    cores = {}
    c = iter(colours)
    c.next()
    richness_data = {'total':data}
    print 'Starting analysis'
    start = time.time()
    co_dets[iteration] = {}
    refinery[iteration] = {}
    zoo[iteration] = {}
    all_collapsed_clusters[iteration] = []

#    colours = ['iz']

    for colour in colours:
        if (colour == colours[-1]) and (single_det == False):
            continue

        #        colour_select = ['iz']
        colour_select = colours

        if colour in colour_select:
            pass
        else:
            print 'Skipping %s'%(colour)
            continue

        start_colour = time.time()
        print '####################################################################################'
        print '                                     %s'%(colour)
        print '####################################################################################'
        try:
            next_colour = c.next()
        except StopIteration:
            print 'This should be the final colour'

        run_function = 'parallel_detect'
        args = [{'source[i].cuts':None},{'source[i].selection':None}]


#        source[colour],job_server = sel_edit2(core.selection,colour+'_width',width[colour],perfect_clusters=False,pp=False,sourcedata=data)
        source[colour] = deepcopy(core)


#        raise SystemExit('@94')
        if iteration in [2,3]:
            source[colour],job_server = sel_edit2(source[colour].selection,colour+'_gradient_function',gradient[colour],perfect_clusters=False,pp=False,sourcedata=data)

        print '%s source selection function created'%(colour)

        # Determine the colour range ORCA will be exploring in C-M space
        intercepts = arange(range[colour][0],range[colour][1],range[colour][2])


#        _intercepts = [intercepts[0],intercepts[-6]]
#        intercepts = _intercepts
#        pp = True
        
        #ATLAS testing for multi-survey flag
#        intercepts = [intercepts[1]]
#        source[colour].selection['colours']['gr']['intercept'] = 1.49147
        
        #g-r = 1.59 = intercept[28] for snapshot testing....
        #intercepts = intercepts[26:30]
#        intercepts = [intercepts[28]]







#        intercepts = intercepts[18:24]
#        intercepts = intercepts[18:20]
#        intercepts = [intercepts[22]]

#        raise SystemExit(0)
#        intercepts = [intercepts[-2]]

        if magmod > 0:
            # This sets an override is we requre the "magmod mode" - this is a progressive scan brightwards in the cmr filter selection, with repeat detections in each.
            # This magmod variable defines how many magnitudes brighter than the limit are explored (ie maglim-magmod) and with what sampling interval (dmag)
            # Will not activate if magmode <= 0

            source[colour].selection['overrides']['depth_mod'] = magmod
            source[colour].selection['overrides']['dmag'] = dmag


        source[colour].selection['defaults']['colour'] = colour
        source[colour].selection['limits']['cmr_range'] = range
        

        
        gladders = False
        if 'gladders' in fileIO.keys() and fileIO['gladders']:
            gladders = True

        if iteration == 1 and do_train:
            # The training cluster is used to define what a red-sequence looks like - it will set (unless subsequently overwritten) the filter slope and width
            # This is done via a CMR fit to a cluster in an ORCA detection object

            print 'Defining training cluster'
            from dataIO import load
            path = train_path+'A2631_evo.cluster'
            # This "evo" uses a spline-fit to define the relationship between cm20 (sequence normalisation) and sequence slope - ie the ORCA cmr filter varies as a function of normalisation


            t_index = 0

            if is_mock:
                path = train_path+'/mock_training.cluster'
                t_index = 0     ## there is often more than one cluster in the detection class: use the 1st one
                

            if os.path.isfile(train_path+'%s.cluster'%(dataset)):
                print 'Using dataset-specific training cluster'
                path = train_path+'%s.cluster'%(dataset)
                
            source[colour].selection['io']['training_cluster'] = 'None'
                
#Colour  gr	 cm20 1.491471	 m -0.048221	 w 0.152384
#Colour  ri	 cm20 0.535871	 m -0.016589	 w 0.152384
#Colour  iz	 cm20 0.332336	 m -0.022785	 w 0.152384

            t_indices = {}
            for _col in colours:
                t_indices[_col] = 0
            try:
                trainer = load(path)
                tclusters = trainer.clusters
                source[colour].selection['io']['training_cluster'] = path
            except:
                trainer = None
                tclusters = []
                print 'Warning: training stipulated, but training set not found!'
            if len(tclusters) > 1:
                from operator import indexOf
                bases = [_cl.basis for _cl in trainer.clusters]
                if not False in [_c in bases for _c in colours]:
                    for _c in colours:
                        t_indices[_c] = indexOf(bases,_c)
                print 'Assigned trainer cluster for each colour'
            failed = []
            for _col in colours:
                t_index = t_indices[_col]
#                if 1:
                try:
                    # This is the attempt to fit the slope, width and intercept (the latter not really used) from the training cluster data
                    source[colour].selection['colours'][_col]['slope'],source[colour].selection['colours'][_col]['width'],source[colour].selection['colours'][_col]['intercept'] = trainer.clusters[t_index].gladders_fit(fit=True,colour=_col,verbose=False,legacy=True,old_method=True)
                    source[colour].selection['colours'][_col]['trainer_width'] = deepcopy(source[colour].selection['colours'][_col]['width'])
                    if fitted_width:
                        source[colour].selection['colours'][_col]['fitted_width'] = source[colour].selection['colours'][_col]['width']
#                else:
                except:
                    failed.append(_col)

            if len(failed) > 0 and (len(colours) != len(failed)):
                seq = deepcopy(colours)
                seq.reverse()
                seq = iter(seq)
                while 1:
                    reddest = seq.next()
                    if reddest in failed:
                        pass
                    else:
                        break

                for _col in failed:
                    source[colour].selection['colours'][_col]['slope'],source[colour].selection['colours'][_col]['width'],source[colour].selection['colours'][_col]['intercept'] = trainer.clusters[t_index].gladders_fit(fit=True,colour=reddest,verbose=False,legacy=True)
                    if fitted_width:
                        source[colour].selection['colours'][_col]['fitted_width'] = source[colour].selection['colours'][_col]['width']

                
            for _col in source[colour].selection['colours'].keys():
                if _col in trainer.selection['colours'].keys() and len(trainer.selection['colours'][_col]['fitfuncs'].keys()) > 0 and do_evo:
                    source[colour].selection['colours'][_col]['fitfuncs'] = trainer.selection['colours'][_col]['fitfuncs']
                    print 'Set learning data for %s'%(_col)


        if source[colour].selection['overrides'].has_key('slope'):
            for key in source[colour].selection['overrides']['slope'].keys():
                source[colour].selection['colours'][key]['slope'] = source[colour].selection['overrides']['slope'][key]

        if source[colour].selection['overrides'].has_key('width'):
            for key in source[colour].selection['overrides']['width'].keys():
                source[colour].selection['colours'][key]['width'] = source[colour].selection['overrides']['width'][key]


        if singular_width:
            maxwidth = max([source[colour].selection['colours'][_col]['width'] for _col in colours])
            print 'All filter widths will be set to the maxwidth [singular_width=True]'

            
            for col in colours:
                source[colour].selection['colours'][col]['width'] = width_factor*maxwidth
                try:
                    source[colour].selection['colours'][col]['fitted_width'] *= width_factor
                except:
                    pass
        else:
            for col in colours:
#                print col,source[colour].selection['colours'][col]['width'], width_factor
                source[colour].selection['colours'][col]['width'] *= width_factor
                try:
                    source[colour].selection['colours'][col]['fitted_width'] *= width_factor
                except:
                    pass
        fw = ''
        if fitted_width:
            fw = '(widths will be auto-adjusted)'

        if do_train:
            print '\nLearnt from training cluster: %s'%(fw)
        else:
            print '\nSet by user: %s'%(fw)
        for col in colours:
            print 'Colour  %s\t cm20 %f\t m %f\t w %f'%(col,source[colour].selection['colours'][col]['intercept'],source[colour].selection['colours'][col]['slope'],source[colour].selection['colours'][col]['width']/width_factor)
        print '\n\n'
##        raise SystemExit(0)

        for col in source[colour].selection['colours'].keys():
            if col != colour:
                source[colour].selection['colours'][col]['width'] = 100.0
          

        cores[colour] = deepcopy(source[colour])
        if colour == colours[-1] and single_det:
            try:
                source[colour].selection['overrides']['return_detection'] = True
            except:
                source[colour].selection['overrides'] = {}
                source[colour].selection['overrides']['return_detection'] = True
            if is_mock and perfect:
                source[colour].selection['overrides']['make_perfect'] = True


#        core.selection['colours']['gr']['width'] = 10.2
#        source[colour].selection['colours']['gr']['width'] = 10.2

#        core.selection['colours']['gr']['width'] = 10.2
#        source[colour].selection['colours']['gr']['width'] = 1.2; intercepts = list(intercepts)+list(intercepts)


#        pp = True
#        intercepts = [1.35]; pp = False
#        intercepts = [1.31,1.35,1.39]; pp = True

        if source[colour].selection['overrides'].has_key('parent_density'):
            from detections import parent_mean_density
            mean_densities,job_server = parent_mean_density(source[colour],intercepts,colour+'_intercept',grid=grid,job_server=job_server,ncpus=ncpus,pp=pp)
            
        if len(intercepts) == 0:
            print 'ERROR: <<intercepts>> has zero-length: no CMR to explore!'
            raise SystemExit(0)
        source[colour],job_server = sel_edit2(source[colour].selection,colour+'_intercept',numpy.array(intercepts),perfect_clusters=False,pp=pp,ncpus=ncpus,sourcedata=data,server=job_server,verbose=True)
        if len(intercepts) == 1:
            source[colour].selection['colours'][colour]['intercept_index'] = 0
        else:
            for i in xrange(len(intercepts)):
                source[colour][i].selection['colours'][colour]['intercept_index'] = i
        print '%s range of detection classes created'%(colour)
        profile.update('%s detection classes created'%(colour))
        if type(source[colour]) != type([]):
            source[colour] = [source[colour]]

        if source[colour][0].selection['overrides'].has_key('parent_density'):
            for i in xrange(len(source[colour])):
                source[colour][i].selection['overrides']['mean_density'] = mean_densities[i]
            print '%s mean density override applied'%(colour)


        completed = []
        from tools import mask

        if pp:
            completed,job_server = pp_run(source[colour],run_function,args,dependent_functions,modules,ncpus=ncpus,idle_delay=1200.0,verbose=True,server=job_server) # Behzad changed idle_delat from 360.0 to 1200.0
        else:
            for s in source[colour]:
                completed.append(parallel_detect(s.cuts,s.selection,plot=False))
        print '%s 1st stage detections completed'%(colour)
        profile.update('%s 1st stage detections completed'%(colour))
        
#        profile.profile()
#        raise SystemExit(0)


        if colour == colours[-1] and single_det:
            source[colour] = completed
        else:
            for slice in xrange(len(completed)):
                source[colour][slice].addclusters(completed[slice],perfect=(is_mock and perfect) and (colour == colours[-1] and single_det))
        print '%s clusters added to detections class'%(colour)

        #from dataIO import save
        #source[colour][0].save('./snapshot_data/snapshot_first_stage.dat')
        #raise SystemExit(0)

        
        if colour != colours[-1] or (single_det == False):
            co_dets[iteration][colour] = accelerate(source[colour],job_server=job_server,perfect=(is_mock and perfect),colour=colour,full_data=data,richness=richness,ncpus=ncpus)
            if richness:
                colour_uber = None
                co_dets[iteration][colour],colour_uber = co_dets[iteration][colour]
                for key in colour_uber.src_data.keys():
                    if not key == 'total':
                        if key in richness_data.keys():
                            print 'Key %d already in richness_data - how?!'
                        richness_data[key] = colour_uber.src_data[key]
                
        else:
            if richness:
                for d in source[colour]:
                    richness_data[d.colour_dna()] = d.cuts
            co_dets[iteration][colour] = source[colour]
            
        profile.update('%s co-detections completed'%(colour))
        print '%s co-detections for all slices completed'%(colour)

#        if pp and (cosma == False):
        if pp:
            cosma = True          #### pp_kill seems to cause problems with some machines (geryon-test)
            if cosma == True:
                try:
                    job_server.destroy()
                except:
                    pass
            else:
                try:
                    pp_kill(job_server)
                except:
                    try:
                        job_server.destroy()
                    except:
                        pass
            del job_server
            job_server = None
            profile.update('%s job_server purged'%(colour))


    print 'Collapsing each colour down to a single detection (splitting for perfect clusters as required)'
    from detections import detections as det
    all_slices = []
    for c in colours:
        try:
            all_slices = all_slices + co_dets[iteration][c]
        except KeyError:
            pass
    perfect_slices = []
    if (is_mock and perfect):
        for sf in all_slices:
            pc_det = det(None,sf.selection,dummy=True)
            pc_det.addclusters(sf.perfect_clusters,verbose=False)
            sf.perfect_clusters = []
            sf.perfect_cluster_galaxies = []
            try:
                pc_det.selection['overrides']['as_perfect'] = True
            except:
                pc_det.selection['overrides'] = {}
                pc_det.selection['overrides']['as_perfect'] = True
            perfect_slices.append(pc_det)
            
    collapsed_clusters = cherrypick(all_slices,sci,metric=core.selection['detection']['metric'],include_solo=True,perfect=False,verbose=True,status_update=True)
    if len(perfect_slices) > 0:
        perfect_clusters = cherrypick(perfect_slices,sci,metric=core.selection['detection']['metric'],include_solo=True,perfect=False,verbose=True,status_update=True)
    else:
        perfect_clusters = False

    profile.update('Collapsed each colour down to single detection')

    base_core = deepcopy(cores[colours[0]])
    base_core.addclusters(collapsed_clusters,perfect=perfect_clusters)
    source = base_core
    if post_proc:
        print 'Post processing clusters for BCG selection'
        base_core.src_data = richness_data
        base_core.post_proc(rfrac=0.5)
        profile.update('Post-processed clusters to locate BCG')
        del base_core.src_data

    if richness:
        print 'Calculating cluster richnesses'
        base_core.src_data = richness_data
        base_core.pp_richness(ncpus=ncpus,uber=False,rfrac=0.5,cosma=cosma)
        profile.update('Cluster richnesses calculated')
        del base_core.src_data

    if fileIO.has_key('compare'):
        print 'Adding external data'
        for key in fileIO['compare'].keys():
            data = fileIO['compare'][key]
            base_core.compare(data['source'],colour=data['colour'],shape=data['shape'],type=data['type'],points=data['points'],name=key)
        profile.update('Added external data')

    if fileIO['save'] != False:
        file = fileIO['save']

        strings = []
        for colour in colours:
            if (colour == colours[-1]) and (single_det == False):
                continue

            from dataIO import save
            prefix = fileIO['prefix']
            field = True
            if fileIO.has_key('zoo') and cluster_zoo:
                save_zoo = fileIO['zoo']
                if save_zoo:
                    zoo_filename = prefix+str(iteration)+'_'+colour+'_zoo.dat'
                    print 'Saving co_detection data to %s'%(zoo_filename)
                    try:
                        save(zoo[iteration][colour],zoo_filename)
                    except:
                        pass
            if fileIO.has_key('co_dets'):
                add_co_dets = fileIO['co_dets']
            if fileIO.has_key('field'):
                field = fileIO['field']
            if field == False:
                base_core.field_galaxies = []
            if add_co_dets:
                cd_filename = prefix+'_'+colour+'_codet.dat'
                print 'Saving co_detection data to %s'%(cd_filename)
                if field == False:
                    for d in co_dets[iteration][colour]:
                        del d.field_galaxies
                        d.field_galaxies = []
                save(co_dets[iteration][colour],cd_filename)

        
        print 'Saving this data to %s'%(file)
        base_core.wrap_toggle(mode='observed')
        save(base_core,file)
        #now remove lock file if it exists....
        break_lock(fileIO)

        base_core.wrap_toggle(mode='detected')
        if fileIO.has_key('symlink') and fileIO['symlink']:
            try:
                os.symlink(file,file.split('/')[-1])
            except:
                pass
    profile.update('Files saved')

    learn = False
    if iteration <= 2 and learn:
        for colour in colours:
            if (colour == colours[-1]) and (single_det == False):
                continue

            if base_core.selection['io'].has_key('observed_data') == False:
                mets = sci.metrics(base_core)
                filtered = sci.filter(base_core,'quality',0.5,1.1)
            else:
                filtered = base_core[colour]
            if (iteration == 1) and len(filtered.clusters) > 0:
                gradient[colour] = filtered.slope_evol(plot=False,colour=colour,perfect=False,passback=False)[0]
            if (iteration == 2) and len(filtered.clusters) > 0:
                width[colour] = filtered.magwidth(plot=False,colour=colour,poly=2,dmag=0.2)

    profile.update('Algorithmic learning')
    stop = time.time()
    return base_core


def status_update(grid,dataset,sleep=5.0):
    global export,cluster_zoo,profile,msgs,subset_index,to_disk,loc,version_id,codets_disk,cluster_zoo,gladders_fit,metric,range,fileIO,stripe,done,init,filecount
    from subprocess import Popen,PIPE

    ss_id = subset_index
    if 1:

        if 1:
            patches = {}
            from tools import define_grid
            import fits_table as fits
            if not 'stripe' in dataset:
                grd,ra_rng,dec_rng = define_grid(105,-0.5,267,68,4,0.75,plot=False)
#                grid = load('/export/disk/dmurphy/data/dr7full/dr7full_grid.dat')
                dr = ra_rng[1]-ra_rng[0]
                cell_factor = 1.0
                fig = plt.figure()
                ax = plt.subplot(111,aspect='equal')
                ncells = len(grid.keys())

            else:
#                grid = load('/export/disk/dmurphy/data/dr7stripes/dr7stripes_grid.dat')
                ra_rng = [360,0]
                dec_rng = [-30,90]                
                dr = 0
                cell_factor = 0.5
                ncells = sum([len(grid[k].keys()) for k in grid.keys()])
                fig = plt.figure(figsize=plt.figaspect(6.1125/15.9))
                ax = plt.subplot(111,aspect='equal')
                plt.subplots_adjust(left=0.05,right=0.98)


            out = ''
            if 1:
                done = []
                import time
                import gobject
                import gtk
                import matplotlib
                matplotlib.use('GTKAgg')
                canvas = ax.figure.canvas

                if  'stripe' in dataset:
                    cr_grid = load('/export/disk/dmurphy/data/dr7full/dr7full_grid.dat')
                    for k in cr_grid.keys():
                        _g = cr_grid[k]
                        rec = plt.Rectangle((_g['xmin'],_g['ymin']),_g['r'],_g['r'],lw=1.0,ec='k',fc='None',ls='dotted',alpha=0.4)
                        ff = plt.gca().add_patch(rec) 
                        


                for k in grid.keys():
                    if 'stripe' in dataset:
                        stripe = k
                        for j in grid[stripe].keys():
                            sid = grid[stripe][j]['sub_id']
                            fn = loc+version_id+'_stripe%dregion%d.dat'%(stripe,sid)
                            if os.path.isfile(fn):
                                done.append(sid)
                                
                            _g = grid[stripe][j]
                            mod  = 0
                            if _g['xc'] > 360.0:
                                mod = -360.0
                            rec = plt.Rectangle((_g['xc']-0.5*cell_factor*_g['r']+mod,_g['yc']-0.5*cell_factor*_g['r']),cell_factor*_g['r'],cell_factor*_g['r'],lw=1.0,ec='k',fc='Red',ls='dotted',alpha=0.4)
                            data = '%s\n%s'%(stripe,_g['sub_id'])
                            plt.text(_g['xc']+mod,_g['yc'],data,horizontalalignment='center',fontsize=6.0,verticalalignment='center')
                            patches[j] = plt.gca().add_patch(rec)

                    else:
                        _g = grid[k]
                        sid = grid[k]['sub_id']
                        fn = loc+version_id+'_region%d.dat'%(sid)
                        if os.path.isfile(fn):
                            done.append(sid)
                        rec = plt.Rectangle((_g['xmin'],_g['ymin']),cell_factor*_g['r'],cell_factor*_g['r'],lw=1.0,ec='k',fc='Red',ls='dotted',alpha=0.4)
                        plt.text(_g['xc'],_g['yc'],str(_g['sub_id']),horizontalalignment='center',fontsize=8.0,verticalalignment='center')
                        patches[k] = plt.gca().add_patch(rec)

                manager = plt.get_current_fig_manager()
                manager.canvas.draw()
                init = True
                tstart = time.time()
                percent = 0.0
                done = []
                from tools import key_lookup,stripe_lookup
                filecount = 0

                def do_update(*args):
                    global done,init
                    if init == False:
                        time.sleep(sleep)
#                    print plt.gcf().get_size_inches()
                    init = False
                    scanner = []
                    logged = []
                    exceptions = []
                    allfiles = os.listdir(loc)
                    allfiles = [loc+f for f in allfiles]
                    if len(allfiles) != filecount:
                        for k in grid.keys():
                            if 'stripe' in dataset:
                                stripe = k
                                for j in grid[stripe].keys():
                                    sid = grid[stripe][j]['sub_id']
                                    fn = loc+version_id+'_stripe%dregion%d.dat'%(stripe,sid)
                                    lf = '%slogs/%s_stripe%dregion_%d.log'%(loc,version_id,stripe,sid)
                                    if fn in allfiles:
                                        scanner.append(j)
                                    elif os.path.isfile(lf):
                                        error_string = 'grep Traceback %s'%(lf)
                                        _exceptions = Popen(error_string.split(),stdout=PIPE).communicate()[0].splitlines()
                                        if len(_exceptions) > 0:
                                            exceptions.append(sid)
                                        else:
                                            logged.append(sid)
                                        

                            else:
                                _g = grid[k]
                                sid = grid[k]['sub_id']
                                fn = loc+version_id+'_region%d.dat'%(sid)
                                lf = '%slogs/%s_region_%d.log'%(loc,version_id,sid)
                                if os.path.isfile(fn):
                                    scanner.append(sid)
                                elif os.path.isfile(lf):
                                    error_string = 'grep Traceback %s'%(lf)
                                    _exceptions = Popen(error_string.split(),stdout=PIPE).communicate()[0].splitlines()
                                    if len(_exceptions) > 0:
                                        exceptions.append(sid)
                                    else:
                                        logged.append(sid)
                                

                    if len(exceptions) > 0:
                        new = exceptions
                        for id in new:
                            stripe = stripe_lookup(grid,id)
                            if not stripe:
                                id = key_lookup(grid,id)
                                _g = grid[id]
                                
                            else:
                                _g = grid[stripe][id]

                            mod  = 0
                            if _g['xc'] > 360.0:
                                mod = -360.0
                            patches[id].set_facecolor('Blue')

                    if len(logged) > 0:
                        new = logged
                        for id in new:
                            stripe = stripe_lookup(grid,id)
                            if not stripe:
                                id = key_lookup(grid,id)
                                _g = grid[id]
                                
                            else:
                                _g = grid[stripe][id]

                            mod  = 0
                            if _g['xc'] > 360.0:
                                mod = -360.0
                            patches[id].set_facecolor('Orange')

                    if len(scanner) > len(done):
                        new = []
                        for id in scanner:
                            if not id in done:
                                new.append(id)
                        for id in new:
                            stripe = stripe_lookup(grid,id)
                            if not stripe:
                                id = key_lookup(grid,id)
                                _g = grid[id]
                                
                            else:
                                _g = grid[stripe][id]

                            mod  = 0
                            if _g['xc'] > 360.0:
                                mod = -360.0
                            patches[id].set_facecolor('Green')

                        done = scanner




                    percent = 100.0*(float(len(done)))/float(ncells)
                    duration = time.time()-tstart
                    if duration > 3600:
                        dh = int(duration /(60.0*60.0))
                        dm = int(60*(duration /(60.0*60.0)-int(dh)))
                        ds = int(duration-((3600.0*dh)+(60.0*dm)))
                        duration = '%dhr:%dm:%ds'%(dh,dm,ds)
                    elif duration > 60.0:
                        dm = duration /(60.0)
                        ds = int(60*(dm-int(dm)))
                        duration = '%dm:%ds'%(dm,ds)
                    else:
                        duration = '%ds'%(int(duration))
                        
                    plt.title('%1.1f percent complete, monitoring for %s'%(percent,duration))
                    

                    manager.canvas.draw()


                    return True
                plt.xlim(max(ra_rng)+dr,min(ra_rng)-dr)
                plt.ylim(min(dec_rng)-dr,max(dec_rng)+dr)
                print plt.xlim()
                print plt.ylim()
                gobject.idle_add(do_update,grid,done)
                plt.xlabel('RA')
                plt.ylabel('DEC')
                plt.show()
                raise SystemExit(0)






def caveats(dataset,core):
    global export,cluster_zoo,profile,msgs,subset_index,to_disk,loc,version_id,codets_disk,cluster_zoo,gladders_fit,metric,range,fileIO,stripe,cg_thresh,pe_thresh
    fileIO = {}
    from tools import which_limits
    from numpy import degrees,radians
    from dataIO import Progress
    import time
    from dataIO import save


    def learning(core,width_learn=False):
        #out = det.learn(plot=True,plot_bins=True,mode='cm20-slope',minbin=[5,5,5,5],s=[1,1,1,1],weight=True);plt.show() ############
        a1 = numpy.array([ 0.74387715,  0.74387715,  0.74387715,  0.74387715,  1.80056787,1.80056787,  1.80056787,  1.80056787])
        a2 = numpy.array([-0.03873526, -0.04950758, -0.04656256, -0.04720442,  0.        ,0.        ,  0.        ,  0.        ])
        a3 = 3
        core.selection['colours']['gr']['fitfuncs']['cm20-slope'] = [a1,a2,a3]
        a1 = numpy.array([ 0.24445163,  0.24445163,  0.24445163,  0.24445163,  1.3142992 ,1.3142992 ,  1.3142992 ,  1.3142992 ])
        a2 = numpy.array([-0.02717082, -0.00116641, -0.02491031, -0.04442966,  0.        ,0.        ,  0.        ,  0.        ])
        a3 = 3
        core.selection['colours']['ri']['fitfuncs']['cm20-slope'] = [a1,a2,a3]
        a1 = numpy.array([ 0.18375215,  0.18375215,  0.18375215,  0.18375215,  0.66862139,0.66862139,  0.66862139,  0.66862139])
        a2 = numpy.array([-0.02593966, -0.0261061 ,  0.00612762, -0.04328547,  0.        ,0.        ,  0.        ,  0.        ])
        a3 = 3
        core.selection['colours']['iz']['fitfuncs']['cm20-slope'] = [a1,a2,a3]
        
        if width_learn:
            a1 = numpy.array([ 0.86900515,  0.86900515,  0.86900515,  0.86900515,  1.74789547,1.74789547,  1.74789547,  1.74789547])
            a2 = numpy.array([ 0.05723322,  0.06872452,  0.06763398,  0.07093645,  0.        ,0.        ,  0.        ,  0.        ])
            a3 = 3
            core.selection['colours']['gr']['fitfuncs']['cm20-scatter'] = [a1,a2,a3]
            a1 = numpy.array([ 0.40715096,  0.40715096,  0.40715096,  0.40715096,  1.07478876,1.07478876,  1.07478876,  1.07478876])
            a2 = numpy.array([ 0.05442641,  0.07610284,  0.05520114,  0.07456444,  0.        ,0.        ,  0.        ,  0.        ])
            a3 = 3
            core.selection['colours']['ri']['fitfuncs']['cm20-scatter'] = [a1,a2,a3]
            

    ########### GLOBAL ################
    cluster_store = export+'%s/'%(dataset)
    if not os.path.isdir(cluster_store):
        print 'Creating cluster store for this dataset:\n\t%s'%(cluster_store)
        if cluster_store[:2] == './':                                                 ## a crappy fix for the python version asts3 is currently running. Looks like this is not needed
            cluster_store = cluster_store.replace('./',os.path.realpath('./')+'/')    ## in Python 2.6.7
        os.makedirs(cluster_store)                                              
    core.selection['detection']['metric'] = metric
    core.selection['limits']['cmr_range'] = range
    st_rep = ''
    smart_rename = True
    sub_sampled = False
    fileIO = {'save':to_disk,'prefix':loc+'%s'%(version_id),'field':False,'co_dets':codets_disk,'zoo':cluster_zoo,'gladders':gladders_fit,'symlink':True}        

    ## fix for annoying Dropbox symlink problem
    if 'Dropbox' in os.getcwd():
        print 'Possibly running in Dropbox directory: disabling symlinks'
        fileIO['symlink'] = False


    
    ########### GLOBAL OVERRIDES ######
    if core.selection.has_key('overrides') == False:
#    if 1:
        core.selection['overrides'] = {}
    core.selection['overrides']['cg_thresh'] = cg_thresh
    core.selection['overrides']['pe_thresh'] = pe_thresh

    ###################################

    if core.selection['detection'].has_key('empty'):
        core.selection['limits'] = {}
        core.clusters = []
        core.cluster_galaxies = []


    ######## OVERRIDE DESCRIPTIONS ########
    ######## THESE ARE NON-STANDARD OVERRIDES THAT EASILY PROPAGATE THROUGH ORCA (AND PP) 
    ######## WITHOUT THE NEED FOR PASSING ARGS TO FUNCTIONS. THEY ARE ADDED TO THE .selection
    ######## CONFIGURATION DICT THAT LIVES WITH EACH DETECTION

    override_info = """\
    keep_matches : when merging clusters in science.cherrypick(), when registering associates to a chosen cluster,
                   keep_matches=True allows the actual associate cluster objects to remain attached to the best cluster
    qhull_bounded: if activated, it will populate random points (via rbox) around the survey to prevent excessive cell 
                   areas at the boundary

    """





    if 1:
        if 's82' in dataset and 'blk' not in dataset and  'mini' not in dataset:
            print '\n\t####################################################################################'
            print '\t####################################################################################'
            print '\t                #now detecting in region %d'%(subset_index)
            print '\t####################################################################################'
            print '\t####################################################################################'

            core.selection['limits']['ymin'],core.selection['limits']['ymax'] = -0.022060094634751632, 0.02206095178515935
            core.selection['limits']['xmin'],core.selection['limits']['xmax'] = which_limits(subset_index)
            core.selection['detection']['subset_index'] = subset_index
            st_rep = '%s_region%d.dat'%(version_id,subset_index)
            if subset_index == -1:
                st_rep = str(version_id)
            else:
                sub_sampled = True
            if sub_sampled:
                fileIO['symlink'] = False

            core.selection['overrides']['keep_matches'] = False
            
        elif 'dr7full' in dataset:
            print '\n\t####################################################################################'
            print '\t####################################################################################'
            print '\t                #now detecting in region %d'%(subset_index)
            print '\t####################################################################################'
            print '\t####################################################################################'
            if core.selection.has_key('grid'):
                core.set_limits()
            width_learn = False
            learning(core,width_learn=False)

            st_rep = '%s_region%d.dat'%(version_id,subset_index)

            #### temp
            sub_sampled = True
#            sub_sampled = False
            ####
            if sub_sampled:
                fileIO['symlink'] = False
                smart_rename = False

            if fitted_width:
                core.selection['overrides']['fitted_width'] = True

                core.selection['overrides']['fitted_width'] = True
                core.selection['overrides']['fitted_width'] = True


#            core.selection['overrides']['area'] = 0.002684466378799872
#            core.selection['overrides']['area'] = 0.0053114007166266677
#            core.selection['overrides']['area'] = 0.0032597091742569873

        elif 'dr7stripes' in dataset:
            print '\n\t####################################################################################'
            print '\t####################################################################################'
            print '\t                #now detecting in stripe %d, sub_region %d'%(stripe,subset_index)
            print '\t####################################################################################'
            print '\t####################################################################################'

            if core.selection.has_key('grid'):
                core.set_limits()
            width_learn = False
            try:
                learning(core,width_learn=False)
            except KeyError:
                pass          ##most likely an empty cell => minimal selection data in core
                
            st_rep = '%s_stripe%dregion%d.dat'%(version_id,stripe,subset_index)

            #### temp
            sub_sampled = True
#            sub_sampled = False
            ####
            if sub_sampled:
                fileIO['symlink'] = False
                smart_rename = False

            if fitted_width:
                core.selection['overrides']['fitted_width'] = True

#            core.selection['overrides']['detection_area'] = 30.92190561804853




        ######## could still be a "generic" sub-sampled dataset (e.g. Stripe 82 data from CFHTLS)
        elif subset_index != -1:
            print '\n\t####################################################################################'
            print '\t####################################################################################'
            print '\t                #now detecting in region %d'%(subset_index)
            print '\t####################################################################################'
            print '\t####################################################################################'

            core.selection['detection']['subset_index'] = subset_index
            st_rep = '%s_region%d.dat'%(version_id,subset_index)
            if subset_index == -1:
                st_rep = str(version_id)
            else:
                sub_sampled = True
            if sub_sampled:
                fileIO['symlink'] = False

            core.selection['overrides']['keep_matches'] = False
            learning(core)

            



        elif 'sas' in dataset:
            if core.selection.has_key('grid'):
                core.set_limits()
                if 'sdss' in dataset:
                    learning(core,width_learn=False)

            core.selection['overrides']['negative_slope'] = ['iz']

            fileIO = {'save':to_disk,'prefix':loc+'%s'%(version_id),'field':False,'co_dets':codets_disk,'zoo':cluster_zoo,'gladders':gladders_fit}


        elif dataset == 'md08':
            if core.selection.has_key('grid'):
                core.set_limits()
            learning(core,width_learn=False)
            core.selection['overrides']['fitted_width'] = True
            core.selection['overrides']['negative_slope'] = ['iy']

            fileIO = {'save':to_disk,'prefix':loc+'%s'%(version_id),'field':False,'co_dets':codets_disk,'zoo':cluster_zoo,'gladders':gladders_fit}


        elif 'cs82blk1' in dataset and subset_index == 1:
            print 'REMOVE THIS CAVEAT - IT WAS FOR TESTING ONLY!!!!!'
            time.sleep(5)
            core.selection['overrides']['area'] = 0.00214517625806


        elif 'usno' in dataset:
            core.selection['colours']['br']['slope'] = -0.075186066


        elif 'dr7' in dataset:
            if maxBCG_mode == False:
                core.selection['overrides']['no_maxBCG'] = True


        elif dataset == 'md08':
            if core.selection.has_key('grid'):
                core.set_limits()
            learning(core,width_learn=False)
            core.selection['overrides']['fitted_width'] = True
            core.selection['overrides']['negative_slope'] = ['iy']

        elif dataset == 'blue82':
            from math import pi
            core.selection['overrides']['area'] = pi*pi*4*(270.0/(360.0*360.0))

        elif dataset == 'herschel':
            core.selection['overrides']['area'] = radians(5.5)*radians(4.0)
            core.selection['overrides']['probthresh'] = 0.024

        elif dataset == 'md05':
            core.selection['overrides']['area'] = 0.0040266995681998083   #v1
            core.selection['overrides']['area'] = 0.0035313516054450702   #v2
        elif dataset == 'md08v2':
            core.selection['overrides']['area'] = 0.0035792885050664968
            #            core.selection['overrides']['probthresh'] = 0.024


        elif dataset == 'd7_str82mini':
            core.selection['overrides']['area'] = 0.01150485590914231




        elif dataset == 'mockv3':
            from math import pi
            batch = True
            export = '/data/rw1/dmurphy/'
            cluster_store = export+'%s/'%(dataset)
            if not os.path.isdir(cluster_store):
                print 'Creating cluster store for this dataset:\n\t%s'%(cluster_store)
                os.makedirs(cluster_store)
            fileIO = {'save':to_disk,'prefix':loc+'%s'%(version_id),'field':False,'co_dets':codets_disk,'zoo':cluster_zoo,'gladders':gladders_fit,'symlink':True}

#            core.selection['overrides']['cg_thresh'] = 0.5
#            core.selection['overrides']['pe_thresh'] = 0.8

            core.selection['overrides']['cg_thresh'] = cg_thresh
            core.selection['overrides']['pe_thresh'] = pe_thresh
            print 'Set cg_thresh, pe_thresh to %f,%f'%(core.selection['overrides']['cg_thresh'],core.selection['overrides']['pe_thresh'])

#            print fileIO['save']
#            raise SystemExit(0)


            learnfile = '%s.learn'%(dataset)
            if os.path.isfile(learnfile):
                learn = load(learnfile)
                for c in core.selection['colours'].keys():
                    if c in learn.selection['colours'].keys() and len(learn.selection['colours'][c]['fitfuncs'].keys()) > 0:
                        core.selection['colours'][c]['fitfuncs'] = learn.selection['colours'][c]['fitfuncs']
                        print 'Set learning data for %s'%(c)
                del learn


        elif dataset == 'mock':
            from math import pi
            xmin,xmax,ymin,ymax = core.selection['limits']['xmin'],core.selection['limits']['xmax'],core.selection['limits']['ymin'],core.selection['limits']['ymax']
            rmax = max([xmin,xmax,ymin,ymax])
            core.selection['overrides']['area'] = pi*rmax*rmax
            core.selection['overrides']['detection_area'] = pi*degrees(rmax)*degrees(rmax)

        elif 'sas' in dataset:
            if core.selection.has_key('grid'):
                core.set_limits()
            if 'sdss' in dataset:
                learning(core,width_learn=False)

            core.selection['overrides']['negative_slope'] = ['iz']

            fileIO = {'save':to_disk,'prefix':loc+'%s'%(version_id),'field':False,'co_dets':codets_disk,'zoo':cluster_zoo,'gladders':gladders_fit}

        elif dataset == 'sa22dxs':
            core.selection['overrides']['negative_slope'] = ['JK']
            core.selection['overrides']['slope'] = {'JK':-0.07}

#            core.selection['overrides']['voronoi area'] = True


        elif dataset == 'stripemini':
            core.selection['overrides']['negative_slope'] = ['gr','ri','iz']

        elif dataset in ['region1','region2']:
            core.selection['overrides']['negative_slope'] = ['gr','ri','iz']



        elif dataset == 'sa22dxsv1':
            core.selection['overrides']['negative_slope'] = ['JK']
            core.selection['overrides']['slope'] = {'JK':-0.07}
#            core.selection['colours']['iJ']['slope'] = -0.022784595709920016            #using i-z sdss slope
#            core.selection['colours']['JK']['slope'] = -0.07

#            core.selection['colours']['iJ']['width'] = 0.152
#            core.selection['colours']['JK']['width'] = 0.152

        elif dataset == 'dxs':
            print 'DXS data - no presets used yet'
#            core.selection['overrides']['voronoi area'] = True

#        if dataset == 's82':
        if dataset == 'eq82':

            core.selection['overrides']['negative_slope'] = ['iz']


            a1 = numpy.array([ 0.86466796,  0.86466796,  0.86466796,  0.86466796,  1.74318448,1.74318448,  1.74318448,  1.74318448])
            a2 = numpy.array([ 0.12726208,  0.18252447,  0.25178933,  0.39670103,0.        ,        0.        ,  0.        ,  0.        ])
            a3 = 3
            core.selection['colours']['gr']['fitfuncs']['cm20-z'] = [a1,a2,a3]

            a1 = numpy.array([ 0.42076046,  0.42076046,  0.42076046,  0.42076046,  1.07764264,1.07764264,  1.07764264,  1.07764264])
            a2 = numpy.array([ 0.17124683,  0.45968845,  0.42745739,  0.54799091,0.        ,        0.        ,  0.        ,  0.        ])
            a3 = 3
            core.selection['colours']['ri']['fitfuncs']['cm20-z'] = [a1,a2,a3]

            a1 = numpy.array([ 0.28580838,  0.28580838,  0.28580838,  0.28580838,  0.54471581,0.54471581,  0.54471581,  0.54471581])
            a2 = numpy.array([ 0.2074385 ,  0.41155726,  0.47409369,  0.56000051,0.        ,0.        ,  0.        ,  0.        ])
            a3 = 3
            core.selection['colours']['iz']['fitfuncs']['cm20-z'] = [a1,a2,a3]

            a1 = numpy.array([ 0.04140306,  0.04140306,  0.04140306,  0.04140306,  0.62359011, 0.62359011,  0.62359011,  0.62359011])
            a2 = numpy.array([-0.04013978, -0.08469228, -0.07210202, -0.04185243,  0.        ,        0.        ,  0.        ,  0.        ])
            a3 = 3
            core.selection['colours']['gr']['fitfuncs']['z-slope'] = [a1,a2,a3]

            a1 = numpy.array([ 0.04140306,  0.04140306,  0.04140306,  0.04140306,  0.62359011,        0.62359011,  0.62359011,  0.62359011])
            a2 = numpy.array([-0.01886636, -0.01885962, -0.03336839, -0.0619371 ,  0.        ,        0.        ,  0.        ,  0.        ])
            a3 = 3
            core.selection['colours']['ri']['fitfuncs']['z-slope'] = [a1,a2,a3]

            a1 = numpy.array([ 0.04140306,  0.04140306,  0.04140306,  0.04140306,  0.62359011,        0.62359011,  0.62359011,  0.62359011])
            a2 = numpy.array([-0.02741922,  0.00092257, -0.02515618, -0.02303591,  0.        ,        0.        ,  0.        ,  0.        ])
            a3 = 3
            core.selection['colours']['iz']['fitfuncs']['z-slope'] = [a1,a2,a3]


        if fileIO['save'] != False:
            ######## if this is a batched submission, might be an identical process also trying to allocate a file.
            ######## for this reason, sleep for a random period between 0-2 secs
          

            repos = cluster_store
#            if sub_sampled:
#                repos = loc
            file = repos+'%s.dat'%(str(version_id))
            if st_rep != '':
                file = repos+st_rep
            lockfile = file.split('.dat')[0]+'.lock'

            if (os.path.isfile(file) or os.path.isfile(lockfile)) and not smart_rename:
                if no_clobber:
                    if not os.path.isfile(file):
                        print '#main(): Filename aready reserved by another ORCA session, and will not clobber - exiting'
                    else:
                        print '#main(): File aready exists, and will not clobber - exiting'
                    raise SystemExit(0)
                else:
                    print '\n!!!!!!!!!Will overwrite file %s'%(file) 
                    print '\t Wating 10 seconds for an abort... (CTRL+\ is safest here)'
                    timer = Progress(10)
                    for i in arange(10):
                        time.sleep(1)
                        timer.update(0)

            elif os.path.isfile(file) and smart_rename:
                if no_clobber:
                    ver = 2
                    while(os.path.isfile(file) or os.path.isfile(lockfile)):
                        root = None
                        if '/' in file:
                            root = '/'.join(file.split('/')[:-1])+'/'
                            file = file.split('/')[-1]
                        splitter = '.dat'
                        if '_v' in file:
                            splitter = '_v'
                            ver = file.split('_v')[1].split('.')[0]
                            ver = int(ver)
                            ver = ver +1
                        file = file.split(splitter)[0]+'_v%d'%(ver)+'.dat'
                        if root != None:
                            file = root+file
                        lockfile = file.split('.dat')[0]+'.lock'    
                        tsleep = 0.5*numpy.random.random()
#                        print '[file:%s] sleeping for %fs'%(file,tsleep)
                        time.sleep(tsleep)
                    if not os.path.isfile(lockfile):
                        os.system('touch %s'%(lockfile))
                        print '#main(): Set lockfile %s'%(lockfile)

                else:
                    print '\n!!!!!!!!!Will overwrite file %s'%(file) 
                    print '\t Wating 10 seconds for an abort... (CTRL+\ is safest here)'
                    timer = Progress(10)
                    for i in arange(10):
                        time.sleep(1)
                        timer.update(0)
            fileIO['save'] = file
            if not os.path.isfile(lockfile):
                os.system('touch %s'%(lockfile))
                print '#main() Set lockfile %s'%(lockfile)


        if dataset == 'mockv3':
#            print fileIO['save']
#            raise SystemExit(0)
            if batch:
                try:
                    mock_lookup = load('mockv3lookup.dat')
                    k = mock_lookup.keys()
                except:
                    mock_lookup = {}

                from dataIO import save
                for key in mock_lookup.keys():
                    if (abs(mock_lookup[key]['cg_thresh']-cg_thresh) < 1e-4) and (abs(mock_lookup[key]['pe_thresh']-pe_thresh) < 1e-4) and os.path.isfile(key):
                        print 'There is already a very similar combination for detection @ %s, bailing'%(key)
                        raise SystemExit(0)
                mock_lookup[fileIO['save']] = {'cg_thresh':core.selection['overrides']['cg_thresh'],'pe_thresh':core.selection['overrides']['pe_thresh']}
                print fileIO['save']
                print mock_lookup[fileIO['save']]
                print cg_thresh,core.selection['overrides']['cg_thresh']
                print pe_thresh,core.selection['overrides']['pe_thresh']
                
#                raise getOut('abort')
                save(mock_lookup,'/data/rw1/dmurphy/PanSTARRS/mockv3lookup.dat')
                print 'Saved mock lookup data'

            

        if is_mock == False:
            fileIO['compare'] = {}
            fileIO['compare']['maxBCG'] = {'source':'./ext_cats/maxbcg.cat','colour':'r','shape':'circle','type':'cluster','points':False}
            fileIO['compare']['Abell'] = {'source':'./ext_cats/acoz.cat','colour':'y','shape':'square','type':'cluster','points':False}
            fileIO['compare']['RASS'] = {'source':'./ext_cats/rass.cat','colour':'b','shape':'square','type':'cluster','points':False}


    if 1:
        profile_plot = False
        info = """\ DR7 stripe runs, with final algorithm parameters

        """
        if subset_index == -1:
            profile_plot = True
            profile_plot = False

        profiler_loc = export+'%s/%s.profile'%(dataset,dataset)
        profiler_path = []
        if not os.path.isfile(profiler_loc):
            profiler = {}
        else:
            profiler = load(profiler_loc)
        if profiler.has_key(version_id) and sub_sampled:
            profiler_path.append(version_id)
            profiler_path.append(subset_index)
            if type(profiler[version_id]) != type({}):
                profiler[version_id] = {}
                save(profiler,profiler_loc)
#        elif profiler.has_key(version_id):
#            new_key = fileIO['save'].split('/')[-1]
        elif sub_sampled:
            profiler_path.append(version_id)
            profiler_path.append(subset_index)
            profiler[version_id] = {}
            save(profiler,profiler_loc)
        else:
            profiler_path.append(version_id)
            profiler = {}
            profiler[version_id] = {}
            save(profiler,profiler_loc)

        del profiler




        fileIO['profile'] = {'plot':profile_plot,'output':True,'info':info,'profiler':profiler_loc,'profile_path':profiler_path}
        print '#main(): caveats routine complete'
                             


    ####################################################################################
    ################################# INPUT ARGUMENTS ##################################
    ####################################################################################

def break_lock(fileData):
    """
    \break_lock: remove the lockfile from the system


    """
    try:
        fname = fileData['save'].split('.dat')[0]+'.lock'
        import os
        cmd = 'rm %s > /dev/null'%(fname)
        os.system(cmd)
    except:
        pass

if __name__ =='__main__':
    global range,intercept,colours,gradient,ncpu,dc,width,co_dets,init_width,fileIO,sci,cores,all_collapsed_clusters,refinery,job_server,source,export,zoo,profile,msgs,subset_index,to_disk,loc,version_id,codets_disk,cluster_zoo,gladders_fit,metric,range,stripe,do_train,cosma,magmod,dmag,cg_thresh,pe_thresh,train_path,do_evo,is_mock,log
    import inspect

    import sys
#    if __debug__:
#        print 'ORCA running under the Python debugger...'

    if len(sys.argv) == 0:
        string = """
    ORCA launcher - loads cluster-finding algorithm (via detect_rapid.py)

    usage:
    =====

    ./detect_rapid.py -dataset [s82] -logged [True] -cat [False] -start [0] -stop [42] -sub_ids x,y,z -reverse [False] -stripe [-1]
                     ...... -stripes i,j,k -version_id [detection] -width_factor [1.0] -ncpu [16] -loc [/export/disk/dmurphy/cat/]
                     ...... -metric [reduced_flux] -nmin [5] -simulate [False]

    switches / parameters:
    =====================
    
    We divide these into the main aspects of ORCA cluster detection

    Survey
    ------
        dataset       : name of the dataset - check detect_generic.py for latest list
        vid           : the "identifier" to use for output products - eg. -vid orca_test_run
        subset_index  : if the data are gridded (see ./grid_data) - run ORCA on this subset_index
        stripe        : if this is SDSS data, define the stripe number (eg: -stripe 82) [default: -1]
        uber          : if activated, then composite "uber" clusters comprising members + associate members will be generated [default: False]
        ncpu          : set the number of CPUs (cores!) to use with ORCA [default: 8, off:1]
    
    Data I/O
    --------
        save          : if activated, will save output to disk
        make_fits     : [redundant]
        save_codets   : [for debugging] save each co-detection file during ORCA detection (will take more disk space)
        save_field    : [redundant]
        cluster_zoo   : [redundant]
        no_clobber    : if filename exists - quit
        getcore       : just write core detection file to disk and then quit [default: False]


    General parameters
    ------------------
        cg_thresh     : what is the merging threshold for systems sharing a fraction cg_thresh of common galaxies?
        pe_thresh     : what is the merging threshold for systems overlapping in projection by a fraction pe_thresh?
        metric        : what defines the "best" cluster for multiple detections? [default: 'reduced flux']
        ngalmin       : what is the minimum number of galaxies in a cluster? [default: 5]
        post_proc     : should we choose a more appropriate BCG from associated cluster members? [default: True]
        richness      : calculate cluster richness? (uses BGC for the moment...) [detault: True]

    Photometric filter
    ------------------
        dc            : the increment in colour when scanning through the CMR [default: 0.04]
        gladders_fit  : apply a fit to the training-cluster sequence as per Gladders++ (1998)
        width_factor  : the multiple by which to multiply the sequence width [default: 0.5]
        maxBCG_mode   : [redundant]
        singular_width: should all filters be set to the same (maximum) width?
        fitted_width  : adopt the width as measured from the training clusters? [default: False]
        single_det    : should the final colour be run on it's own?
        trainer       : use trainer cluster for filter parameters
        evoln         : (if available from trainer) use sequence learning?
        magmod        : number of mags brighter than limit to probe (-ve: disables)
        dmag          : mag increment with which to peel back faintwards

    Voronoi parameter
    ------------------
        rho_crit      : the minimum density threshold before the FoF-stage of cluster growth halts
        probthresh    : the probability threshold (Kiang 1966)
        d_pt          : [debug] set the increment range for Probthresh for parameter exploration 
        d_rc          : [debug] set the increment range for rho_crit for parameter exploration 


    Mock parameters
    ---------------
        mock          : Is the source data a mock catalogue?    
        perfect       : Also create "perfect" clusters during run? (ie, from DM-halo data)
        halomass      : filter all halos below a specified halomass [default: (10**) 13.0]
        outlier_metric: define what quantity is used to filter out any non-member outliers
        z             : filter all halos above specified redshift [default: 0.0]
        dz            : 



    eg.
    ==

    ./detect_rapid.py -dataset herschel -metric flux
    ./detect_rapid.py -dataset snapshot -vid 'snapshot_tests' -ncpu 8 -trainer True -magmod -1 -mock False

        """
        os.system('clear')
        print string

        raise SystemExit(0)


    t1 = time.time()


    argv,dataset = rdarg(argv,"-dataset",str,'orca_test')#
#    argv,dataset = rdarg(argv,"-dataset",str,'cfhtls_v2')
#    argv,dataset = rdarg(argv,"-dataset",str,'snapshot')#
#    argv,dataset = rdarg(argv,"-dataset",str,'test_100sd_aper3')#
#    argv,dataset = rdarg(argv,"-dataset",str,'a226')#
#    argv,dataset = rdarg(argv,"-dataset",str,'ATLAS_SGC6_A5')#
#    argv,dataset = rdarg(argv,"-dataset",str,'ATLAS_test_aper5')#
#    argv,dataset = rdarg(argv,"-dataset",str,'mock')
#    argv,dataset = rdarg(argv,"-dataset",str,'cfhtls_orcatest')
#    argv,dataset = rdarg(argv,"-dataset",str,'mock_conito')

    argv,subset_index = rdarg(argv,"-subset_index",int,-1)#
    argv,stripe = rdarg(argv,"-stripe",int,-1)#
    argv,ncpu = rdarg(argv,"-ncpu",int,8)#
    

    argv,version_id = rdarg(argv,"-vid",str,'orca_test_run')#
#    argv,version_id = rdarg(argv,"-vid",str,'ATLAS_SGC6_A5_run')#
#    argv,version_id = rdarg(argv,"-vid",str,'base_run')#
#    argv,version_id = rdarg(argv,"-vid",str,'atlas_tests')#
    argv,uber_mode = rdarg(argv,"-uber",bool,True)#

    # data I/O
    argv,to_disk = rdarg(argv,"-save",bool,True)#
    argv,to_fits = rdarg(argv,"-make_fits",bool,False)#
    argv,codets_disk = rdarg(argv,"-save_codets",bool,False)#
    argv,field_disk = rdarg(argv,"-save_field",bool,False)#
    argv,cluster_zoo = rdarg(argv,"-cluster_zoo",bool,False)#
    argv,no_clobber = rdarg(argv,"-no_clobber",bool,True)#
    argv,getcore = rdarg(argv,"-getcore",bool,False)#

    # photometric filter parameters
    argv,dc = rdarg(argv,"-dc",float,0.04)#
    argv,gladders_fit = rdarg(argv,"-gladders_fit",bool,True)#
    argv,width_factor = rdarg(argv,"-width_factor",float,0.5)#
    argv,maxBCG_mode = rdarg(argv,"-maxBCG_mode",bool,False)#
    argv,singular_width = rdarg(argv,"-singular_width",bool,False)#
    argv,fitted_width = rdarg(argv,"-fitted_width",bool,False)    
    argv,single_det = rdarg(argv,"-single_det",bool,True)#
    argv,do_train = rdarg(argv,"-trainer",bool,1)#
    argv,do_evo = rdarg(argv,"-evoln",bool,1)#
    argv,magmod = rdarg(argv,"-magmod",float,-1)#
    argv,dmag = rdarg(argv,"-dmag",float,0.2)#


    # voronoi parameters
    argv,rho_crit = rdarg(argv,"-rho_crit",float,2.5)#
    argv,probthresh = rdarg(argv,"-probthresh",float,0.0125)#
    argv,d_pt = rdarg(argv,"-d_pt",float,0.0)#
    argv,d_rc = rdarg(argv,"-d_rc",float,0.0)#

    #mock detector parameters
    argv,is_mock = rdarg(argv,"-mock",bool,False)#
    argv,perfect = rdarg(argv,"-perfect",bool,True)#
    argv,halomass = rdarg(argv,"-halomass",float,13.0)#
    argv,outlier_metric = rdarg(argv,"-outlier_metric",str,'radius')#
    argv,outlier_no_clobber = rdarg(argv,"-outlier_no_clobber",bool,False)   ## ?
    argv,redshift = rdarg(argv,"-z",float,0.0)#
    argv,dredshift = rdarg(argv,"-dz",float,1.9)#
    
    #general detector parameters

    argv,cg_thresh = rdarg(argv,"-cg_thresh",float,0.2)#
    argv,pe_thresh = rdarg(argv,"-pe_thresh",float,0.2)#
    argv,metric = rdarg(argv,"-metric",str,'reduced_flux')#
    argv,smallest_cluster = rdarg(argv,"-ngalmin",int,5)#
    argv,post_proc = rdarg(argv,"-post_proc",bool,True)#
    argv,richness = rdarg(argv,"-richness",bool,False)#



#    log_loc = getcfgProp('orca','logging',cfgloc='./',addslash=False,verbose=False)
    
    from sys import argv
    print 'Running this process on %s. Started on %s'%(socket.gethostname(),time.ctime())


    ## Do some critical checks before proceeding:
    from tools import whereis   
    from subprocess import Popen,PIPE
    qhull_cmd = getcfgProp('orca','qhull',cfgloc='./',addslash=False,verbose=False)
    rbox_cmd = getcfgProp('orca','rbox',cfgloc='./',addslash=False,verbose=False)
    search_string = 'compute convex hulls and related structures'
    print qhull_cmd
    if qhull_cmd and (not os.path.isfile(qhull_cmd)):
        print 'The orca.cfg defined qhull does not exist, looking elsewhere...'
        qhull_cmd = None
    if qhull_cmd:
        qhull_cmd = [qhull_cmd]
        rbox_cmd = [rbox_cmd]
        print qhull_cmd
        p1 = Popen(qhull_cmd, stdout=PIPE,stderr=PIPE)
        output = p1.communicate()[0]
        if not 'compute convex hulls and related structures' in output:
            print 'Default qhull does not work, looking in binary repository...'
            qhull_cmd[0] = None

    else:
        print 'No qhull referenced in orca.cfg, searching paths:'
        qhull_cmd = whereis('qhull')
        qhull_cmd = [qhull_cmd]
        rbox_cmd = whereis('rbox')
        rbox_cmd = [rbox_cmd]
        print qhull_cmd
        if qhull_cmd[0] != None:
#            qhull_cmd = [qhull_cmd]
            print qhull_cmd
            p1 = Popen(qhull_cmd, stdout=PIPE,stderr=PIPE)
            output = p1.communicate()[0]
            if not 'compute convex hulls and related structures' in output:
                print 'Path-located qhull does not work, looking in binary repository...'
                qhull_cmd[0] = None


    if qhull_cmd[0] == None:
        print 'Searching binary repository:'
        bins = os.listdir('./lib/')
        qhulls = []
        rboxes = []
        for b in bins:
            if 'qhull' in b:
                print 'Located %s as qhull binary'%('./lib/%s'%(b))
                qhulls.append('./lib/%s'%(b))
                rboxes.append('./lib/%s'%(b.replace('qhull','rbox')))
        if len(qhulls) == 0:
            print 'No alternative qhull binaries found: please install'
            print 'the qhull / rbox paths can be added to orca.cfg:'
            print 'qhull                : /home/binaries/qhull'
            raise SystemExit(0)
        found = False
        for qhull_cmd in qhulls:
            qhull_cmd = [qhull_cmd]
            p1 = Popen(qhull_cmd, stdout=PIPE,stderr=PIPE)
            output = p1.communicate()[0]
            if 'compute convex hulls and related structures' in output:
                found = True
                break
            
        if found:
            qhull_cmd = qhull_cmd[0]
            print 'Using %s'%(qhull_cmd)
            setcfgProp('orca','qhull',qhull_cmd,cfgloc='./')
            setcfgProp('orca','rbox',qhull_cmd.replace('qhull','rbox'),cfgloc='./')

        elif qhull_cmd == None:
            print 'No viable qhull binaries found: please install'
            print 'the qhull / rbox paths can be added to orca.cfg:'
            print 'qhull                : /home/binaries/qhull'
            raise SystemExit(0)

        


        else:
            print 'None of the alternative qhull binaries work: please install'
            print 'the qhull / rbox paths can be added to orca.cfg:'
            print 'qhull                : /home/binaries/qhull'
            raise SystemExit(0)

        

    else:
        #the whereis qhull works - let's just stick it in the .cfg file...
        if not getcfgProp('orca','qhull',cfgloc='./',addslash=False,verbose=False):
            if type(qhull_cmd) == type([]):
                qhull_cmd = qhull_cmd[0]
                rbox_cmd = qhull_cmd.replace('qhull','rbox')
            setcfgProp('orca','qhull',qhull_cmd,cfgloc='./')
            setcfgProp('orca','rbox',rbox_cmd,cfgloc='./')
            

#    raise SystemExit(0)
    






#<<<<<<< TREE #behzad 31OCT2017
    #    argv,dataset = rdarg(argv,"-dataset",str,'cfhtls_v2')
#    argv,dataset = rdarg(argv,"-dataset",str,'snapshot')#
#    argv,dataset = rdarg(argv,"-dataset",str,'test_100sd_aper3')#
#    argv,dataset = rdarg(argv,"-dataset",str,'a226')#
    argv,dataset = rdarg(argv,"-dataset",str,'orca_test')#

#    argv,dataset = rdarg(argv,"-dataset",str,'ATLAS_test_aper5')#

#    argv,dataset = rdarg(argv,"-dataset",str,'mock')
#    argv,dataset = rdarg(argv,"-dataset",str,'cfhtls_orcatest')
#    argv,dataset = rdarg(argv,"-dataset",str,'mock_conito')

    argv,subset_index = rdarg(argv,"-subset_index",int,-1)#
    argv,stripe = rdarg(argv,"-stripe",int,-1)#
    argv,ncpu = rdarg(argv,"-ncpu",int,8)#
    

    argv,version_id = rdarg(argv,"-vid",str,'orca_test_run')#
#    argv,version_id = rdarg(argv,"-vid",str,'base_run')#
#    argv,version_id = rdarg(argv,"-vid",str,'atlas_tests')#
    argv,uber_mode = rdarg(argv,"-uber",bool,True)#

    # data I/O
    argv,to_disk = rdarg(argv,"-save",bool,True)#
    argv,to_fits = rdarg(argv,"-make_fits",bool,False)#
    argv,codets_disk = rdarg(argv,"-save_codets",bool,False)#
    argv,field_disk = rdarg(argv,"-save_field",bool,False)#
    argv,cluster_zoo = rdarg(argv,"-cluster_zoo",bool,False)#
    argv,no_clobber = rdarg(argv,"-no_clobber",bool,True)#
    argv,getcore = rdarg(argv,"-getcore",bool,False)#

    # photometric filter parameters
    argv,dc = rdarg(argv,"-dc",float,0.04)#
    argv,gladders_fit = rdarg(argv,"-gladders_fit",bool,True)#
    argv,width_factor = rdarg(argv,"-width_factor",float,0.5)#
    argv,maxBCG_mode = rdarg(argv,"-maxBCG_mode",bool,False)#
    argv,singular_width = rdarg(argv,"-singular_width",bool,False)#
    argv,fitted_width = rdarg(argv,"-fitted_width",bool,False)    
    argv,single_det = rdarg(argv,"-single_det",bool,True)#
    argv,do_train = rdarg(argv,"-trainer",bool,1)#
    argv,do_evo = rdarg(argv,"-evoln",bool,1)#
    argv,magmod = rdarg(argv,"-magmod",float,-1)#
    argv,dmag = rdarg(argv,"-dmag",float,0.2)#


    # voronoi parameters
    argv,rho_crit = rdarg(argv,"-rho_crit",float,2.5)# org_base_run 5  # new_base_run 2.5 # orca_original 10
    argv,probthresh = rdarg(argv,"-probthresh",float,0.0125)# org_base_run 0.025 # new_base run 0.0125 # orca_original 0.01
    argv,d_pt = rdarg(argv,"-d_pt",float,0.0)#
    argv,d_rc = rdarg(argv,"-d_rc",float,0.0)#

    #mock detector parameters
    argv,is_mock = rdarg(argv,"-mock",bool,False)#
    argv,perfect = rdarg(argv,"-perfect",bool,True)#
    argv,halomass = rdarg(argv,"-halomass",float,13.0)#
    argv,outlier_metric = rdarg(argv,"-outlier_metric",str,'radius')#
    argv,outlier_no_clobber = rdarg(argv,"-outlier_no_clobber",bool,False)   ## ?
    argv,redshift = rdarg(argv,"-z",float,0.0)#
    argv,dredshift = rdarg(argv,"-dz",float,1.9)#
    
    #general detector parameters

    argv,cg_thresh = rdarg(argv,"-cg_thresh",float,0.2)# base_run 0.2
    argv,pe_thresh = rdarg(argv,"-pe_thresh",float,0.2)# base_run 0.2
    argv,metric = rdarg(argv,"-metric",str,'reduced_flux')#
    argv,smallest_cluster = rdarg(argv,"-ngalmin",int,5)#
    argv,post_proc = rdarg(argv,"-post_proc",bool,True)#
    argv,richness = rdarg(argv,"-richness",bool,False)#

#=======              #behzad 31OCT2017  
#>>>>>>> MERGE-SOURCE #behzad 31OCT2017 (conflict detected as DAVID updated a different part of the code and his version didn't have my modifications here)
    ####################################################################################
    ############################## SELECTION PARAMETERS ################################
    ####################################################################################

    env_args = {'probthresh':probthresh,'smallest_cluster':smallest_cluster,'redshift':redshift,'dredshift':dredshift,'rho_crit':rho_crit,\
            'halomass':halomass,'d_pt':d_pt,'d_rc':d_pt,'dc':dc,'metric':metric,'outlier_metric':outlier_metric,'outlier_no_clobber':outlier_no_clobber}
    cosma = False
    node = None
    if socket.gethostname() == 'asts3.phyast.dur.ac.uk' or socket.gethostname() == 'icc-graphics.phyast.dur.ac.uk':
        export = '/export/disk/dmurphy/'
        if socket.gethostname() == 'icc-graphics.phyast.dur.ac.uk':
            if ncpu > 1:
                print 'Changing ncpu=8'
                ncpu = 8

    elif socket.gethostname()[0] == 'm':
        #could this be a cosma node?
        host = socket.gethostname()
        try:
            node = int(host[1:])
            cosma = True
            
        except ValueError:
            #seems not...
            pass

    elif socket.gethostname() in ['pollux','castor','m3001','m3002']:

        cosma = True

        if socket.gethostname() in ['pollux','castor','m3001','m3002']:
            cosma = True
        else:
            print 'I am running on cosma4 node %s'%(socket.gethostname())
        export = '/data/rw1/dmurphy/'
        if ncpu > 1:
            ncpu = ncpu

    elif 'dubris' in socket.gethostname():
        export = '/export/home/dmurphy/'
        if ncpu > 1:
            print 'Changing ncpu=8'
            ncpu = 8

    else:
        export = '/export/disk/dmurphy/'


    note = """now disregarding platform-dependent locations for output - will read from orca.cfg
              setting outpath to '/export/disk/dmurphy/' in orca.cfg is tanamount to the above.
               """

    export = getcfgProp('orca','outpath',cfgloc='./',addslash=True,verbose=False)
    train_path = getcfgProp('orca','trainpath',cfgloc='./',addslash=True,verbose=False)
    if train_path == None or os.path.isdir(train_path) == False:
        #prevent ORCA from using a training cluster...
        if os.path.isdir(train_path) == False:
            print 'WARNING: train path not found!!!'

        do_train = False

    ####################################################################################
    ############################## PARALLEL PYTHON INIT ################################
    ####################################################################################
    ncpus = ncpu
    pp = True
    if ncpus == 0 or ncpus == 1: pp = False
    completed = []
    #ppservers=("*",)              # use for autodiscovery of clusters nodes
    ppservers = ()
    sci = science()
    dependent_functions = (deploy_vertices,qhull,parallel_detect,sel_edit,sci.richness)
    modules = ('science','parallel_detect','clusters','numpy','parallel_qhull','parallel_vertices','vperc','time','detections')
    args = [{'source[i].cuts':None},{'source[i].selection':None}]
    run_function = 'parallel_detect'

    ####################################################################################
    ################################ DATA SELECTION ####################################
    ####################################################################################

    myt1 = time.time()
    
    ####################################################################################
    # Staging area for variable overrides - be sure to remove these to reset to defaults
    ####################################################################################


    profile = Progress(0,timeline=True,writesteps=True)
    print '#main(): profiler started'
    msgs = None
#    loc = export+'cat/'
    loc = export

    if 1:
        do_wrap = False
        ss_id = None
        grid = None
        stripeN = None

        if 'eq82' in dataset:
            stripeN = 82
            ss_id = subset_index


        if 's82' in dataset:
            from tools import which_limits
            xmin,xmax = which_limits(subset_index)
            do_wrap = xmin < 0
            ss_id = subset_index

        if 'dr7stripes' in dataset:
            ss_id = subset_index
            from tools import define_grid
            import fits_table as fits
            gridfile = 'mygrid.dat'
            gridfile = '/export/disk/dmurphy/data/dr7stripes/dr7stripes_grid.dat'

            grid = load(gridfile)
            if stripe == -1:
                status_update(grid,dataset)
#                print 'Please enter a stripe number!'
                raise SystemExit(0)
            else:
                try:
                    grid = grid[stripe]
                except KeyError:
                    raise SystemExit('Cannot find stripe %d in grid. Valid stripes are :\n%s'%(stripe,grid.keys()))
            stripeN = stripe

        if 'dr7full' in dataset or ('dr7' in dataset and len(dataset) > len('dr7')) and 'stripe' not in dataset:
            ss_id = subset_index
            from tools import define_grid
            import fits_table as fits
            grid,ra_rng,dec_rng = define_grid(105,-0.5,267,68,4,0.75,plot=False)
            grid = load('/export/disk/dmurphy/data/dr7full/dr7full_grid.dat')

            out = ''
            if ss_id == -1:
                status_update(grid,dataset)
                raise SystemExit(0)
                for id in grid.keys():
                    out = out + '%d,'%(grid[id]['sub_id'])
                out = out[:-1]
                print '\n\nFull set of %d subset_indices:'%(len(grid.keys()))
                print out+'\n\n'
            do_wrap = False

        elif 'dr7' in dataset:
            ss_id = subset_index
            do_wrap = True
            if 'maxBCG' in version_id:
                maxBCG_mode = True
        
        else:
            if subset_index != -1 and subset_index > 0:
                print 'Non-standard subset_index set - please be sure this is what you want to do'
                print '(setting subset_index to %s)'%(subset_index)
                ss_id = subset_index



#        data,junk,junk,range,core = getData('',args=env_args,dataset=dataset,wrap=do_wrap,selection=None,subset_index=ss_id,plot=False
        data,junk,junk,range,core = getData('',args=env_args,dataset=dataset,wrap=do_wrap,selection=None,subset_index=ss_id,plot=False,maxBCG_mode=maxBCG_mode,grid=grid,stripe=stripeN,md5=True)
                     
        print '#main(): Core detection object produced'

        try:
            caveats(dataset,core)
        except:
            import sys
            print fileIO
            break_lock(fileIO)
            exc_type,exc_error,trace = sys.exc_info()
            print 'Something happened    :'
            import traceback
#            print aa
#            print '@ line number: '
            traceback.print_exception(exc_type,exc_error,trace)

            traceback.print_tb(trace)
            raise SystemExit('caught')
          
        if data == None:
            try:
                core.save(fileIO['save'])
            except:
                pass
            raise SystemExit(0)



        sci.range = core.selection['limits']['cmr_range']
        colours = core.selection['detection']['sequence']
#        colours = ['br']
        init_width = 0.05
        width = {}
        for c in colours:
            width[c] = init_width


#        colours = ['gr']
#        colours = ['gr','ri']
#        colours = ['ri','iz']
#        colours = ['iz','iy']


        print core.selection['limits']['xmin'],'.....',core.selection['limits']['xmax']
#        print fileIO['prefix']
        print fileIO['save']

#        print fileIO
#        break_lock(fileIO)
#        raise SystemExit(0)


        save(core,'%s_core.dat'%(str(dataset)))
        if getcore:
            break_lock(fileIO)
            raise SystemExit(0)
        core.cuts = data
        cores = {}
        iterations = [1]
        results = {}
        gradient = {}
        co_dets = {}
        refinery = {}
        zoo = {}
        all_collapsed_clusters = {}
        job_server = None
        gc.collect()
        
        for iteration in iterations:
            t2 = time.time()
            if 1:
#            try:
                results[iteration] = run(data,core,iteration)
            else:
#            except:
                print 'Run failed'
                break_lock(fileIO)
                raise SystemExit(0)
            t3 = time.time()

    all = results[1]
    if to_fits and (dataset != 's82' or subset_index == -1):
        from tools import cat
        job_server = cat(all,output='%s_collapsed.dat'%(version_id),job_server=job_server)
        profile.update('output catalogue produced')

    if fileIO.has_key('profile'):

        profile.info = fileIO['profile']['info']
        key = fileIO['profile']['profile_path']
        profiler = load(fileIO['profile']['profiler'],verbose=False)
        if len(key) > 1:
            key,sub_key = key
            profiler[key][sub_key] = profile
        else:
            key = key[0]
            profiler[key] = profile
        save(profiler,fileIO['profile']['profiler'],verbose=False)
        print 'Profiler updated'
        del profiler


        if fileIO['profile']['plot']:
            profile.profile(plot=True)
            plt.title('Runtime profile for dataset: %s, version: %s'%(dataset,version_id))
        elif fileIO['profile']['output']:
            profile.profile()

    det = all
    print
    print '%d galaxies placed into %d clusters'%(len(det.cluster_galaxies),len(det.clusters))
    try:
        from numpy import median
        zz = [cl.z for cl in det.clusters]
        print 'min/median/max redshift is %1.3f / %1.3f / %1.3f'%(min(zz),median(zz),max(zz))
    except:
        pass
    print 'Fin'
#    det.save('dtest.dat')

try:
    plt.show()
except:
    pass

