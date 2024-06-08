import clusters
import socket
from numpy import arange,take,where,array,zeros,std
#import detections as classes
from operator import countOf,indexOf,contains
from copy import copy,deepcopy
import gc
import fits_table as fits
try:
    from lib.pp import pp as ppserver
except:
    pass
import time
#from tools import *
import thread
import sys
from math import degrees
import os
if socket.gethostname() != 'dubris':
    import matplotlib.pyplot as plt
    import matplotlib


class getDataError(Exception):
    import os

    def __init__(self,msg):
        print '\t#########################################################################'
        print '\t\n\n\n\n\t\t%s\n\n\n\n'%(msg)
        print '\t#########################################################################'
        if not 'dubris' in socket.gethostname() and 'Zuul' in msg:
            try:
                import Image
                dog = Image.open('zuul.jpg')
                dpi = plt.rcParams['figure.dpi']
                figsize = dog.size[0]/dpi, dog.size[1]/dpi
                plt.ax = plt.axes([0,0,1,1], frameon=False)
                plt.ax.set_axis_off()
                plt.imshow(dog, origin='lower')
                plt.show()
            except:
                pass
        elif not 'dubris' in socket.gethostname() and 'need' in msg and 'input' in msg:
            try:
                import Image
                j5 = Image.open('../johnny5.png')
                dpi = plt.rcParams['figure.dpi']
                figsize = j5.size[0]/dpi, j5.size[1]/dpi
                plt.ax = plt.axes([0,0,1,1], frameon=False)
                plt.ax.set_axis_off()
                plt.imshow(j5, origin='lower')
                plt.show()
            except:
                pass

        raise SystemExit(0)


class Progress:
    def __init__(self,total,widget=False,timer=True,timeline=False,writesteps=False,pidcheck=None):
        import time
        if type(total) == type(''):
            #assume this is a request for a PID update...
            pidcheck = total
        
        if pidcheck != None:
            import os
            if os.path.isfile(pidcheck):
                self.restore(pidcheck)
                self.profile()
                return


        self.timeline = timeline
        self.widget = widget
        if self.widget:
            self.widgetlock = thread.allocate_lock()
        self.timer = timer
        self.writesteps = writesteps
        if self.writesteps:
            import os
            self.pid_fname = 'profile.%d'%(os.getpid())


        if self.timer: self.start = time.time()
        if timeline == True:
            self.init = time.ctime().split()
            self.timeline = True
            self.time_sequence = []
            self.time_line = []
            self.end = self.start
            return

        print '\n'
        try:
            self.total = len(total)
        except:
            self.total = total
        self.lock = thread.allocate_lock()
        self.count = -1
        self.percent = 0
        self.update(self.total)


    def restore(self,file):
        from dataIO import load
        input = load(file,convert=False)
        self.time_sequence = input['time_sequence']
        self.time_line = input['time_line']
        self.init = input['init']
        self.end = input['start']
        self.start = input['end']
        self.timeline = input['timeline']

    def profile(self,plot=False,store=False,filename=None,total_only=False):
        if not self.timeline:
            return
        if total_only == False:
            print '\nProfile started at : \t\t\t%s'%(' '.join(self.init))
            print 'Step\t Time(s)\t Time(m:s)\tInfo'
        for i in xrange(len(self.time_sequence)):
            t_secs = self.time_sequence[i]
            mins = int((t_secs)/60.0)
            secs = int(60*(((t_secs)/60.0)-int(((t_secs)/60.0))))
            msg = str(self.time_line[i])
            if total_only == False:
                print '%d\t %.1f\t\t %dm:%ds  \t%s'%(i+1,t_secs,mins,secs,msg)
        t_secs = self.end - self.start
        mins = int((t_secs)/60.0)
        secs = int(60*(((t_secs)/60.0)-int(((t_secs)/60.0))))
        msg = 'Total run time'
        total = t_secs
        if total_only:
            print '%dm:%ds  \t%s\n'%(mins,secs,msg)
        else:
            print '%d\t %.1f\t\t %dm:%ds  \t%s\n'%(i+2,t_secs,mins,secs,msg)

        if plot:
            fig = plt.figure()
            indices = arange(len(self.time_sequence))+1
            labels = [str(i) for i in indices]
            fracs = [100.0*(t/t_secs) for t in self.time_sequence]
            plt.pie(fracs,labels=labels,autopct='%1.2f%%',shadow=True)
            plt.title('Runtime profile')

        if store:
            from dataIO import save
            date = self.init
            t0 = date[3].split(':')
            file = '%s_%s_%s@%s_%s.profile'%(date[2],date[1],date[-1],t0[0],t0[1])
            if filename == None:
                filname = file
            output = {}
            output['time_sequence'] = self.time_sequence
            output['time_line'] = self.time_line
            output['init'] = self.init
            output['start'] = self.end
            output['end'] = self.start
            save(output,filename)

    #the callback function
    def update(self,dummy):
        if self.timeline == True:
            t = time.time()
            if len(self.time_sequence) == 0:
                self.time_sequence.append(t-self.start)
            else:
                tslu = t - self.end

                self.time_sequence.append(tslu)
            self.time_line.append(dummy)
            self.end = t
            if self.writesteps:
                output = {}
                output['time_sequence'] = self.time_sequence
                output['time_line'] = self.time_line
                output['init'] = self.init
                output['start'] = self.end
                output['end'] = self.start
                output['timeline'] = self.timeline
                save(output,self.pid_fname,verbose=False)
            return

        #dummy is the job output....
        # must use lock here because += is not atomic
#        print '\t\n'

        self.lock.acquire()
        try:
            if dummy == -1:
                pass
            else:
                self.count += 1
        except:
            raise IOError()

        self.lock.release()
        percent = int(100.0*float(self.count)/float(self.total))
        self.percent = percent
        finish = 100
        progress = '{'+percent*'#'+(finish-percent)*' '
        progress += '}'
        remaining = self.total-self.count
        if remaining == 0:
            #kill the thinger if it's running
            self.stop = True
        if self.timer:
            mins = int((time.time()-self.start)/60.0)
            secs = int(60*(((time.time()-self.start)/60.0)-int(((time.time()-self.start)/60.0))))
            bar = progress+str(percent)+'%'+' %dm:%ds'%(mins,secs)+' %d remain'%(remaining)
        else : bar = progress+str(percent)+'%'
        print bar,'\r',
        self.bar = bar
        sys.stdout.flush()
        if self.count == self.total: print '\t\n'
        if self.widget and self.count ==0:
            threadlock = thread.allocate_lock()            
            thread.start_new_thread(self.thinger,(threadlock.locked(),))

    def thinger(self,status):
        self.widgetlock.acquire()
        stop = False
        oldput = '|'
        while 1:
            if self.percent >= 99: stop=True
            if oldput == '|': output = '/'
            if oldput == '/': output = '-'
            if oldput == '-': output = '\\'
            if oldput == '\\': output = '|'
            print self.bar[:self.percent+1]+output+self.bar[self.percent+2:],'\r',

            if stop:
                self.widgetlock.release()
                thread.exit_thread()
            sys.stdout.flush()
            time.sleep(0.05)
            oldput = output


def run_cmd2(string,stdout=True):
    import os
    import sys
    print 'Iam running %s [stdout=%s] '%(string,str(stdout))
#    os.system('touch /export/disk2/tiles/stripe82/mytester.dat')
    if stdout == True:
        print 'Doing...'
        os.system(string)
    else:
        os.system(string+' >& /dev/null')
    print 'Done - here is the output...'
    sys.stdout.flush()
    return


def run_cmd(string,stdout=False,verbose=True):
    import os
    import sys
    from subprocess import Popen,PIPE
    stdout = True
    if verbose:
        print 'Iam running %s'%(string)
#    os.system('touch /export/disk2/tiles/stripe82/mytester.dat')
    cmd = string.split()
    p1 = Popen(cmd, stdout=PIPE,stderr=PIPE)
#    header = p1.communicate()[0]
    header = p1.communicate()
    if stdout:
        if verbose:
            print 'Done - here is the output...'
            print header[0]
            print header[1]
            sys.stdout.flush()

    return header



def pp_run(source,function,pp_args,dependent_functions,modules,ncpus=8,ppservers=(),idle_delay=1200.0,background_processes=3,verbose=False,server=None,progress=True,protocol='auto',secret=None):
    # Behzad changed idle_delat from 520.0 to 1200.0
    #got zombies? Kill them and/or parents via:
    #ps -ef | grep defunct | grep [username] | grep -v grep | cut -b9-20 | xargs kill -9


    global threadcount,threadlock,bar,percent,thread,run_widget,done
    import os
    import numpy
    from inspect import getargvalues, stack, currentframe
    function_args = getargvalues(currentframe())[0]
    
    #posname, kwname, myargs = getargvalues(stack()[1][0])[-3:]
    #function_args = myargs['locals'].keys()
    ##check that no function args have the same name

    pp_arg_names = [p.keys()[0] for p in pp_args]
    forbidden = []
    for _arg in pp_arg_names:
        ignores = ['verbose','ncpus']
        if _arg in function_args:
            if (_arg in ignores) == False:
                forbidden.append(_arg)
    if len(forbidden) > 0:
        for _arg in forbidden:
            print '#dataIO.pp_run(): ERROR: supplied argument <%s> cannot be used as it can be contaminated by pp_run function args!'%(_arg)
        raise SystemExit(0)

    
    #    from detect import update
# if the function needs to be imported, the module containing it must be the first module in the modules list.
# does it need to be imported?
#    if verbose:print '\tUsing fancy generic submitter...'

#    if socket.gethostname() in ['pollux','castor','m3001','m3002'] or socket.gethostname()[0] == 'm':

#    if socket.gethostname() in ['pollux','castor','m3001','m3002']:



    if socket.gethostname() in ['xx']:
        srv = {}
        srv['castor'] = "172.17.100.21"
        srv['pollux'] = "172.17.100.22"
        srv['m3001'] = "172.17.100.1"
        srv['m3002'] = "172.17.100.2"
        secret = 'orca'
        ppservers = ("172.17.101.1","172.17.100.21","172.17.101.2")
        ppservers = ("172.17.*",)

        ppservers = ("172.17.101.136",)
        ppservers = ("172.17.101.*",)
        ppservers = ("*",)
        ppservers = ("172.17.101.66",)
        ppservers = ("172.17.101.*",)
        ppservers = tuple(["172.17.101.%d"%(i) for i in arange(3,221)])
    string = 'X.'
    try:
        i = 0
        X = eval(function)
    except NameError:
        try:
            X = __import__(modules[0])
        except ImportError:
            if verbose: print '\tModule %s does not exist'%(modules[0])
        try:
            string = 'X.'
            X = eval(string+function)
        except AttributeError:
            if verbose: print '\t Scanning through modules for function <%s>'%(function)
            for mod in modules:
                try:
                    X = __import__(mod)
                    string = 'X.'
                    X = eval(string+function)
                except:
                    if mod == modules[-1]:
                        if verbose:
                            print '\tCannot find function <%s> in modules list:\n %s'%(function,modules)
                            print '\tReturning gracefully (programme will likely fail from this point on)'
                        return []
                    pass
                else:
                    if verbose: print '\tFound function <%s> in module <%s>'%(function,mod)
                    break

    # need to check that the source list hasn't been nested by mistake (often happens):
    i = 0
#    exploded = []
    key = None
    for j in xrange(len(pp_args)):
        if pp_args[j][pp_args[j].keys()[0]] == None and '[i].' in pp_args[j].keys()[0]:
            key = pp_args[j].keys()[0]

    if key != None:
        source_copy = deepcopy(source)
        try:
            tester = eval(key)
        except AttributeError:
            if type(source[0]) == type([]):
                source = source[0]
                try:
                    tester = eval(key)
                except:
                    print '\tToo much nesting....'
                    return []
                source = []
                for i in xrange(len(source_copy)):
                    source.append(source_copy[i][0])
                print '\t#Adjusted input source format'


    if ncpus < 1: pp = False
    else: pp = True
    done = 0
    total = len(source)

#    source[0][0]['wait'] = True

    dummy = 0
    run = []

    proto = 0
    if protocol == 'auto':
        import pickle
        proto = pickle.HIGHEST_PROTOCOL
    
    if server != None:
        job_server = server
        if job_server.get_ncpus() != ncpus:
            job_server.set_ncpus(ncpus)
            
    else:   
        print ppservers
        print 'secret = %s'%(secret)
        job_server = ppserver.Server(ncpus,ppservers=ppservers,proto=proto,secret=secret)
    if verbose: print '\t#Submitting %d jobs to parallel job server. Using %d cores to run <<%s>>'%(len(source),min(ncpus,len(source)),function)
    if progress:
        status = Progress(len(source),widget=False)
    else:
        status = None
        start = time.time()
#    job = ppserver.Template(job_server, X, dependent_functions, modules,callback=status.update)
    argarray = []
    for i in xrange(len(source)):
        a = []
        for j in xrange(len(pp_args)):
            key =  pp_args[j].keys()[0]
            if pp_args[j][key] == None:
                a.append(eval(key))

            elif '[i]' in key and pp_args[j][key] != None:
                a.append(pp_args[j][key][i])

            elif pp_args[j][key] == 'None':
                a.append(None)

            else: a.append(pp_args[j][key])

        a = tuple(a)
        argarray.append(a)
        

        if status:
            run.append(job_server.submit(X,a,dependent_functions,modules,callback=status.update))
        else:
            run.append(job_server.submit(X,a,dependent_functions,modules))

    timer = 0.
    culprit = None
    
    intra_job = -1.
    intra_job_thresh = idle_delay
    jobs_done = 0
    t0 = time.time()
    re_runs = []
    while countOf([r.finished for r in run],True) != len(run):
        time.sleep(1)
        if status:

            status.update(-1)
            if status.count != jobs_done:
                #reset the timer flags
                jobs_done = status.count
                intra_job = time.time()
                hw = False
                ij_thresh = False
            else:
                if intra_job > 0 and (time.time()-intra_job > 10.0) and (time.time()-intra_job < 11.5):
                    _david = 1
#                    print 'intra-Clock started (10s)'

                if intra_job > 0 and (time.time()-intra_job > 0.5*intra_job_thresh) and (hw == False):
                    print '\nPP idle state for half of idle_delay'
                    hw = True

                if intra_job > 0 and (time.time()-intra_job > intra_job_thresh) and (ij_thresh == False):
                    print '\nIntra-job threshold reached - do something!'
                    #strictly speaking, should kill job here and re-supbmit
                    ij_thresh = True
                    state = numpy.array([r.finished for r in run])
                    fails = numpy.where(state == False)[0]
                    print '\nRun numbers %s are still pending - will try to re-run with new ppserver...'%(fails)
                    #kill current attempts first on original server
                    job_server.destroy()
                    job_server = ppserver.Server(ncpus,ppservers=ppservers,proto=proto,secret=secret)
                    #init job_server2
                    job_server2 = ppserver.Server(ncpus,ppservers=ppservers,proto=proto,secret=secret)
                    #run
                    new_server = False
#                    while countOf([r.finished for r in run],True) != len(run):
#                        state = numpy.array([r.finished for r in run])
#                        fails = numpy.where(state == False)[0]
                        


                    for i in fails:
                        a = argarray[i]
                        if status:
                            re_runs.append(job_server2.submit(X,a,dependent_functions,modules,callback=status.update))
                        else:
                            re_runs.append(job_server2.submit(X,a,dependent_functions,modules))
                    re_runtimer = time.time()
                    hw_warn = False
                    while countOf([r.finished for r in re_runs],True) != len(re_runs):
                        time.sleep(1)
                        if ((time.time()-re_runtimer) > 0.5*intra_job_thresh) and (hw_warn == False):
                            print '\nWarning: pp-retry idle state for half of idle_delay on secondary job server'
                            hw_warn = True
                        if (time.time()-re_runtimer) > intra_job_thresh:
                            #kill the job and throw a tantrum
                            print '\nAttempted pp-retry but jobs still exceed idle_delay on secondary job server!'
                            raise SystemExit(0)
                    print 'pp-reruns finished - insering into runs'
                    for i in xrange(len(fails)):
                        run[fails[i]] = re_runs[i]
                    #kill job_server2 to free memory (see below)

        scrap = """\
        if countOf([r.finished for r in run],True)+1 == len(run) and timer == 0 and len(run) > 1:
            timer = time.time()
#            print '\t\nWaiting on final job, starting timer'
        if timer != 0 and time.time()-timer > idle_delay and countOf([r.finished for r in run],True) != len(run):
            print 'Runs now exceeded idle time....'
            job_status = [r.finished for r in run]
            nfinished = countOf(job_status,True)
            if nfinished == len(run):
                break
            print '\tThe remaining task has exceeded idle_delay - throwing exception for pdb to catch'
            workers = [w for w in job_server._Server__workers]
            failed_run = indexOf(job_status,False)
            print '\tDetermined failed_run'
            try:
                culprit = workers[indexOf([w.is_free for w in workers],False)]
                pid = culprit.pid
            except:
                pid = 0

            import os
            print '\tFound pid, imported OS'
            try:
                pp_kill(job_server)
            except IOError:
                pass            
            try:
                job_server.destroy()
                print '\tDestroyed pp server'
            except:
                pass


            path = '/proc/%d'%(pid)
            print '\tWorker process (pid=%d) removed : %s'%(pid,str(os.path.isdir(path)))
            job_server2 = ppserver.Server(ncpus,ppservers=ppservers,proto=proto)
            print '\tProduced new job server'

            if status:
                run[failed_run] = job_server2.submit(X,argarray[failed_run],dependent_functions,modules,callback=status.update)
            else:
                run[failed_run] = job_server2.submit(X,argarray[failed_run],dependent_functions,modules)
            print '\tWaiting for failed job to complete'
            job_server2.wait()
            job_server = job_server2

            
    if culprit == None:
        job_server.wait()
#        try:
#            pp_kill(job_server)
#        except IOError:
#            pass            
#        try:
#            job_server.destroy()
#        except:
#            pass

    else:
        pass
#        try:
#            pp_kill(job_server2)
#        except IOError:
#            pass
#            
#        try:
#            job_server2.destroy()
#        except:
#            pass
    """    
    completed = []
    for job in run:
        job()
        completed.append(job())

    if len(re_runs) > 0:
        print 'Removing secondary job_server used for re-runs'
        try:
            pp_kill(job_server2)
        except IOError:
            print 'pp_kill did not work on job_server destruction'
            
        try:
            job_server2.destroy()
        except:
            print '.destroy() did not work on job_server destruction'
            pass
        del job_server2

    if progress:
        if verbose: print '\t#This job took a total of %d seconds'%(int(time.time()-status.start))
    else:
        if verbose: print '\t#This job took a total of %d seconds'%(int(time.time()-start))
    return completed,job_server

def clean_table():
    from subprocess import Popen,PIPE
    cmd1 = 'ps -l'
    cmd2 = 'head -1'
    cmd1 = cmd1.split()
    cmd2 = cmd2.split()
    p1 = Popen(cmd1, stdout=PIPE)
    p2 = Popen(cmd2, stdin=p1.stdout, stdout=PIPE)
    header = p2.communicate()[0].split()
    pid = indexOf(header,'PID')
    parent_id = indexOf(header,'PPID')

    p1 = Popen(cmd1, stdout=PIPE)
    all_processes = p1.communicate()[0].splitlines()
    all_processes = all_processes[1:]
    for proc in all_processes:
        if 'defunct' or ('wait' and 'sh') in proc:
            proc = proc.split()
    


def pp_kill(job_server,verbose=False):
    """\
    pp_kill(pp.Server) -> void
    
    \t Kills all worker processes in the pp.Server instance.
    This clears the broken pipes that stack up in /proc/<pid>/fd
    when creating many pp.Server instances.

    (New 11/07/09)
    Added signalling control to remove the zombie processes that
    pile up when the workers are closed in this way. There's
    probably a better way to do this.

    Will likely be redundant after pp v1.5.7"""

#    import signal
#    import os
    from subprocess import Popen,PIPE
#    handler = signal.getsignal(signal.SIGCHLD)
#    signal.signal(signal.SIGCHLD, signal.SIG_IGN)
#    import time
#    thispid = os.getpid()

    children = []
    for worker in job_server._Server__workers:
        children.append(worker.pid)

    cmd1 = 'ps -l'
    cmd2 = 'head -1'
    cmd1 = cmd1.split()
    cmd2 = cmd2.split()
    p1 = Popen(cmd1, stdout=PIPE)
    p2 = Popen(cmd2, stdin=p1.stdout, stdout=PIPE)
    header = p2.communicate()[0].split()
    pid = indexOf(header,'PID')
    parent_id = indexOf(header,'PPID')
    parents = []

    for child in children:
        cmd1 = 'ps -l'.split()
        cmd2 = ('grep %d'%(child)).split()
        cmd3 = 'grep python'.split()
        p1 = Popen(cmd1, stdout=PIPE)
        p2 = Popen(cmd2, stdin=p1.stdout, stdout=PIPE)
        p3 = Popen(cmd3, stdin=p2.stdout, stdout=PIPE)
        this_child = p3.communicate()[0].split()
        try:
            parents.append(int(this_child[parent_id]))
        except:
            pass
        p1.wait()
        p2.wait()
        p3.wait()

    import os


#    junk = """\
#    handler = signal.getsignal(signal.SIGCHLD)
#    signal.signal(signal.SIGCHLD, signal.SIG_IGN)
#    import time
    thispid = os.getpid()

    if verbose: os.system('ps -l | head -1')
    for worker in job_server._Server__workers:
        child_pid = worker.pid
        if verbose: os.system('ps -l | grep %d'%(child_pid))
        worker.is_free = False
        if verbose: os.system('ps -l | grep %d'%(child_pid))
        worker.t.send('EXIT')
        if verbose: os.system('ps -l | grep %d'%(child_pid))
        worker.t.close()
        if verbose: os.system('ps -l | grep %d'%(child_pid))

    if verbose: os.system('ps -l')
#    """
    #now just to be sure, try to manually kill all the processes, starting with the parents:

    for parent in parents:
        try:
            os.kill(parent,9)
        except OSError:
#            print '\tProcess %d already dead'%(parent)
            if verbose: print '\tProcess %d already dead'%(parent)
            #already dead
            pass
    for child in children:
        try:
            os.kill(child,9)
        except OSError:
            #already dead
            if verbose: print '\tProcess %d already dead'%(child)
#            print '\tProcess %d already dead'%(child)
            pass
    try:
        if verbose: print '\tWaiting for zombies to die....'
#        os.wait()
    except OSError:
        #this means all the child processes are dead.
        pass

#    signal.signal(signal.SIGCHLD, handler)
    return

def make_core(selection=None,perfect=False,data=None,getdata=False,sim_col=True):
    import detections as classes
    if selection == None:
        selection = {}
        selection['io'] = {}
        selection['io']['sourcefile'] = ''
        selection['io']['dataset'] = ''
        selection['detection'] = {}
        selection['limits'] = {}
        selection['limits']['xmax'] = 1.0
        selection['limits']['xmin'] = 0.0
        selection['limits']['ymax'] = 1.0
        selection['limits']['ymin'] = 0.0

        if sim_col:
            selection['detection']['sequence'] = ['gr']



    observed = selection.has_key('observed_data')
    cuts = {}
    cuts['x'] = []
    cuts['y'] = []
    if data == None and getdata==True:
        data = getData(selection['io']['sourcefile'],obs=observed,dataset=selection['io']['dataset'])
        cuts = slice_data(data,make_mask(data,selection))
    core = classes.detections(cuts,selection,make_perfect=perfect,sourcedata=data)
    return core
    

def slice_data(data,mask):
    xcut = take(data['x'][:],mask)
    ycut = take(data['y'][:],mask)
    try:
        vx = take(data['vx'][:],mask)
        vy = take(data['vy'][:],mask)
        vz = take(data['vz'][:],mask)
        sm = take(data['stellar_mass'][:],mask)
        hm = take(data['halo_mass'][:],mask)
        hid = take(data['halo_id'][:],mask)
        hjm = take(data['halo_jm'][:],mask)
        hjmtree = take(data['halo_jmtree'][:],mask)
        hzcos = take(data['halo_zcos'][:],mask)
        hidrep = take(data['halo_idrep'][:],mask)
    except TypeError:
        vx,vy,vz,sm,hm,hid,hjm,hjmtree,hzcos,hidrep = None,None,None,None,None,None,None,None,None,None

    gcut = take(data['gmag'][:],mask)
    rcut = take(data['rmag'][:],mask)
    icut =  take(data['imag'][:],mask)
    zcut =  take(data['zmag'][:],mask)
    Ycut =  take(data['ymag'][:],mask)

    gr = gcut-rcut
    ri = rcut-icut
    iz = icut-zcut
    zy = zcut-ycut

    cutdata = {'x':xcut,'y':ycut,'stellar_mass':sm,'halo_mass':hm,'halo_id':hid,\
               'halo_jm':hjm,'halo_jmtree':hjmtree,'halo_zcos':hzcos,'gmag':gcut,\
               'rmag':rcut,'halo_idrep':hidrep,'imag':icut,'zmag':zcut,'ymag':Ycut,\
               'vx':vx,'vy':vy,'vz':vz,'gr':gr,'ri':ri,'iz':iz,'zy':zy}
    
    return cutdata


def save(detection,name='output.dat',slimline=True,verbose=True,mode='cPickle',protocol='auto'):
    import numpy
    from pickle import PicklingError
    if mode == 'cPickle':
        import cPickle as pickle
    else:
        import pickle
        
    obj_type = None

    from detections import detections,cluster,galaxy
    dummy = detections(None,None,None,dummy=True)
    dummycl = cluster([None],None,None,dummy=True)
    dummygal = galaxy(-1,None,None,dummy=True)
    js = None

    if type(detection) == type(dummy) and slimline:
        obj_type = 'detection'
        detection.source_data = None
        try:
            js = detection.job_server
            detection.job_server = None
            print 'Stripped job_server for save'
        except:
            pass

    else:
        if type(detection) == type(dummycl):
            obj_type = 'cluster'
        elif type(detection) == type(dummygal):
            obj_type = 'galaxy'

        try:
            detection.parent.source_data = None
            detection.neighbours = None
        except:
            pass

    file = open(name,'w')
    
    proto = pickle.HIGHEST_PROTOCOL
    if protocol != 'auto':
        proto = protocol

    try:
        pickle.dump(detection,file,protocol=proto)
    except:
       
        def pickle_check(obj):
            parent_data = None
            _key = None
            try:
                keys = obj.__dict__.keys()
            except:
                return None
            for key in keys:
                if key == 'parent':
                    parent_data = obj.__dict__[key]
                    del obj.parent
                    continue
                try:
                    daughter = obj.__dict__[key]
                    pick = pickle.dumps(daughter)
#                except:
                except PicklingError:
                    _key = key
                    print '%s failed to pickle'%(_key)
                    if parent_data:
                        obj.parent = parent_data
                        
                    return _key,daughter
            if parent_data:
                obj.parent = parent_data
            return _key

        print 'Picklerror - diagnosing what element must be removed'
        offender = 'start'
        obj = detection
        objs = [detection]
        paths = []
        offenders = []
        while offender != None:
            offender = pickle_check(obj)
            if offender != None:
                offender,obj = offender
                offenders.append(offender)
                if type(obj) in [type([]),type(numpy.array([]))]:
                    tests = numpy.array([pickle_check(_o) != None for _o in obj])
                    if True in tests:
                        failed = numpy.where(tests)[0]
                        print 'items %s failed in %s'%(failed,offender)
                        obj = obj[indexOf(tests,True)]

        print offenders
        if 1:
            return
        obj = detection
        iterator = [obj]
        objs = iter(iterator)
        offender = 'pickle'
        parent_offender = ''
        while offender != None:
            offender = pickle_check(objs.next())
            if offender != None:
                offender,obj = offender
                offenders.append(offender)
                if len(iterator) > 1:
                    print 'Item %d failed in %s'%(indexOf(iterator,obj),parent_offender)
                else:                  
                    iterator = [obj]
                    objs = iter(iterator)

                if type(obj) in [type([]),type(numpy.array([]))]:
                    tests = numpy.array([pickle_check(_o) != None for _o in obj])
                    iterator = obj
                    objs = iter(iterator)
                    parent_offender = offender


                    if True in tests:
                        failed = numpy.where(tests)[0]
                        print 'items %s failed in %s'%(failed,offender)
                        obj = obj[indexOf(tests,True)]

        print offenders





            

    file.close()
    if verbose: print '\t# dataIO : Saved data to %s'%(name)
    if js != None:
        detection.job_server = js


def load(name='output.dat',verbose=True,convert=True,mode='cPickle'):
    if mode == 'cPickle':
        import cPickle as pickle
    else:
        import pickle

    import os
    from clusters import galaxy,cluster
    import numpy
    from detections import detections as det
    dummy = det(None,None,dummy=True)
    dummy_cl = cluster([],-1,dummy=True)
    dummy_gal = galaxy(-1,0,None,None,dummy=True)
#    dummy_cl = cluster([galaxy(-1,0,0,dummy=True)],-1,dummy=True)
    
    if os.path.isfile(name) == False:
        name = '/export/disk/dmurphy/cat/'+name
        if verbose:
            print 'Could not find file, trying alternate location'
    if os.path.isfile(name):
        file = open(name,'r')
        detection = pickle.load(file)
        file.close()
    else:
        if verbose:
            print 'Could not find file %s'%(name)
            return None
    if convert:
        if detection.__class__ == dummy.__class__:
#        if type(detection) == type(dummy):
            # a detection class
            if not detection.selection.has_key('io'):
                detection.convert(verbose=verbose,all=True)
            try:
                detection.wrap_toggle(mode='detected')
            except:
                print 'WARNING - DELIBERATELY SABOTAGED TO PREVENT WRAP (is this OK?)'
                pass




#        elif type(detection) == type(dummy_cl):
        elif detection.__class__ == dummy_cl.__class__:
            # a cluster
            parent = None
            try:
                parent = detection.parent
            except:
                pass

            if parent:
                if not parent.selection.has_key('io'):
                    parent.convert(verbose=verbose,all=True)

            try:
                mag=detection.magnitudes
            except:
                try:
                    detection.convert(verbose=verbose)
                except AttributeError:
                    print 'Not a convertible type (detection, cluster...) - re-loading with convert = False'
                    _dat = load(name,convert=False)
                    return _dat
                    
        elif type(detection) == type([]) or type(detection) == type(numpy.array([])):
            if detection[0].__class__ == dummy.__class__:
#            type(detection[0]) == type(dummy):
                # a list of detection_classes
                for d in detection:
                    if not d.selection.has_key('io'):
                        d.convert(verbose=verbose,all=True)
            
            elif detection[0].__class__ == dummy_cl.__class__:
                # a list of clusters
                parents = []
                for c in detection:
                    parents.append(c.parent)
                parents = unique(parents)
                for d in parents:
                    if d != None:
                        if not d.selection.has_key('io'):
                            d.convert(verbose=verbose,all=True)


    if verbose:
        print '\t# dataIO : loaded data from %s'%(name)
    return detection
    
def make_detections(parameter,colour,cherrypick,index,data,selections,perfect_clusters=True,rlimit=24.0,set_density=False):
    import detections as classes
    observed = False
    fast_forward = [0,1,7]
    if type(data) == type('f'):
        if 'observed_data' in data:
            data = data.split('observed_data')[0]
            observed = True
        data = getData(data,obs=observed,dataset=selections[0]['dataset'],selection=selections[0]) 

    if not selections and type(parameter) != type(False):
        colours = colour
        parameters = parameter
        selections = []
        for k in range(len(colours)):
            a = [colours[j][0] for j in range(len(colours))]
            for m in range(len(a)):
                a[m] = [a[m][n] for n in range(len(a[m]))]

        try:
            z = (parameters[0][5]-parameters[0][6],parameters[0][5]+parameters[0][6],None)
        except:
            z = (None,None,None)

        junk = """\
        print 'Experiment'
        _sel = selections[0]
        print _sel
        print type(_sel)
        """
        if _sel.has_key('custom_box'):
            xmin,xmax = _sel['custom_box']['xlimits']
            ymin,ymax = _sel['custom_box']['ylimits']
            _sel['xmin'] = xmin
            _sel['xmax'] = xmax
            _sel['ymin'] = ymin
            _sel['ymax'] = ymax
#            print 'Overriding limits from    dx=[%f,%f] dy=[%f,%f]'%(parameters[0][2]-parameters[0][4],parameters[0][2]+parameters[0][4],parameters[0][3]-parameters[0][4],parameters[0][3]+parameters[0][4])
            print 'to                        dx=[%f,%f] dy=[%f,%f]'%(xmin,xmax,ymin,ymax)
#            return
        else:
            xmin,xmax = parameters[0][2]-parameters[0][4],parameters[0][2]+parameters[0][4]
            ymin,ymax = parameters[0][3]-parameters[0][4],parameters[0][3]+parameters[0][4]

        print xmin,xmax
        print ymin,ymax
        
        xy = [xmin,xmax,ymin,ymax]
        probthresh = parameters[0][0]
        smallest_cluster = parameters[0][1]
        rho_crit = parameters[0][7]
        a.append(z)
        a.append(xy)
        a.append(rlimit)
        a.append(probthresh)
        a.append(smallest_cluster)
        a.append(rho_crit)
        a.append(data['sourcefile'])
        a.append((degrees(parameters[0][4])*2)**2)
        a.append(parameters[0][8])
        sel = selection_fn(a)
        selections.append(sel)
        parameter = False
   # list of indices for which the full-blown mask-creation can be copied - ie, those
    # parameters having no effect on the mask / cuts (eg, probthresh...)
    cuts = []
    #
    # produce the masks
    #
    _sel = selections[0]
    if _sel.has_key('custom_box'):
        xmin,xmax = _sel['custom_box']['xlimits']
        ymin,ymax = _sel['custom_box']['ylimits']
        _sel['xmin'] = xmin
        _sel['xmax'] = xmax
        _sel['ymin'] = ymin
        _sel['ymax'] = ymax

    if parameter and contains(fast_forward,index):
        print '\t#Creating galaxy cuts'
        c = slice_data(data,make_mask(data,selections[0]))
        cuts.append(c)
        for d in xrange(len(selections)-1):
            cuts.append(copy(c))
    else:
        cuts = [slice_data(data,make_mask(data,selections[i])) for i in xrange(len(selections))]
    if len(cuts) > 1:
        print '\t#Produced %d galaxy cuts'%(len(cuts))
    results = []
    if ((parameter and contains(fast_forward,index)) or ((not parameter and not colour))) and not cherrypick:
        det = classes.detections(cuts[0],selections[0],make_perfect=perfect_clusters,sourcedata=data,set_density=set_density)
        results.append(det)
        for d in xrange(len(cuts)-1):
            det = classes.detections(cuts[d+1],selections[d+1])
            det.perfect_clusters = copy(results[0].perfect_clusters)
            results.append(det)
    else:
        for i in xrange(len(cuts)):
            det = classes.detections(cuts[i],selections[i],make_perfect=perfect_clusters,sourcedata=data,cherrypick=True,set_density=set_density)
            results.append(det)

    if len(results) > 1:
        print '\t#Produced %d detection classes'%(len(results))
    else:
        gr = results[0].selection['gr_intercept']
        ri = results[0].selection['ri_intercept']
        iz = results[0].selection['iz_intercept']
        zy = results[0].selection['zy_intercept']
        ngals = len(results[0].cuts['x'])
        print '\t#Made c(gr,ri,iz,zy)=(%7.5f,%7.5f,%7.5f,%7.5f) detection class with %d galaxies'%(gr,ri,iz,zy,ngals)
    del cuts
    del selections
    return results


def sel_edit2(selection,variable,range,sourcedata=False,cherrypick=False,colour=False,pp=True,rlimit=24.0,perfect_clusters=False,ncpus=8,set_density=False,server=None,verbose=False):
    import fits_table as fits
#    import fits_table as fits
    import numpy
    from detections import detections as det
#    from detections import detections as det
    #1 how many detections to be output?
    core = det(None,None,dummy=True)
    selections = []
    detections = []
    if ((type(range) != type([])) and (type(range) != type(numpy.array([]))) or (len(range) == 1)):
        detections.append(deepcopy(core))
        selections.append(deepcopy(selection))
        if type(range) != type([]) and type(range) != type(numpy.array([])):
            range = [range]
    elif len(range) > 1:
        for r in range:
            detections.append(deepcopy(core))
            selections.append(deepcopy(selection))
    v0,v1,v2 = variable,'',''
    if '_' in variable:
        v0,v1 = variable.split('_')
        if v0 in selection['detection']['sequence']:
            if v1 in selection['colours'][v0].keys():
                for i in xrange(len(detections)):
                    selections[i]['colours'][v0][v1] = range[i]
                    detections[i].selection = selections[i]
                
            else:
                print 'Property %s for colour %s does not exist'%(v1,v0)
                return

        
        elif v0 in selection['limits']['magnitudes'].keys():
            if v1 in selection['limits']['magnitudes'][v0].keys():
                for i in xrange(len(detections)):
                    selections[i]['colours'][v0][v1] = range[i]
                    detections[i].selection = selections[i]
            else:
                print 'Property %s for magnitude %s does not exist'%(v1,v0)
        else:
            print 'Not sure where this goes, putting it into the base location'
            for i in xrange(len(detections)):
                selections[i][variable] = range[i]
                detections[i].selection = selections[i]
            
    else:
        if selection['limits'].has_key(v0):
            for i in xrange(len(detections)):
                selections[i]['limits'][v0][v1] = range[i]
                detections[i].selection = selections[i]
            
        elif selection['detection'].has_key(v0):
            for i in xrange(len(detections)):
                selections[i]['detection'][v0][v1] = range[i]
                detections[i].selection = selections[i]

        elif selection['mock_data'].has_key(v0):
            for i in xrange(len(detections)):
                selections[i]['mock_data'][v0][v1] = range[i]
                detections[i].selection = selections[i]
        
    if len(detections) == 1:
        detections = detections[0]
        detections.get_input(data=sourcedata,perfect=perfect_clusters)

    elif pp:
        args = [{'source[i]':None},{'perfect':perfect_clusters}]
        dets = detections
        import detections as classes
        import detections
        dependent_functions = (selection_fn,slice_data,make_mask,take,make_detections,classes,getData,fits,detections)
        modules = ('clusters','numpy','parallel_qhull','operator','fits_table','operator','detections')
        detections,server = pp_run(dets,'get_input',args,dependent_functions,modules,ncpus=ncpus,idle_delay=1200.0,server=server,verbose=verbose) # Behzad changed idle_delat from 90.0 to 1200.0

    elif not pp:
        for d in detections:
            d.get_input(data=sourcedata,perfect=perfect_clusters)
    
    return detections,server


def sel_edit(selection,variable,range,sourcedata=False,cherrypick=False,colour=False,pp=True,rlimit=24.0,perfect_clusters=True,ncpus=8,set_density=False,server=None):
    from numpy import arange,take,where,array
    import detections as classes
    from operator import countOf,indexOf,contains
    from copy import copy,deepcopy
    import gc
#    import fits_table as fits
    import time
    print '\t# Splitting the selection functions'
    dependent_functions = (selection_fn,slice_data,make_mask,take,make_detections,classes,getData,fits)
    modules = ('clusters','numpy','parallel_qhull','operator','fits_table','operator','detections')
    keywords = ['intercept','slope','width']
    index = -1
    parameter = True
    cherrypick,colour = False,False
    observed = False
    for key in keywords:
        if key in variable:
            cherrypick=True
            colour = True
            parameter = False
    selections = []
    classes = []
    if selection['io'].has_key('observed_data'):
        observed=True

    if not sourcedata:
        if not pp:
            sourcedata = getData(selection['sourcefile'],obs=observed,dataset=selection['dataset'])
    if pp:
        sourcedata = selection['sourcefile']
        if observed:
            sourcedata = sourcedata+'observed_data'


    if type(range) != type(arange(0,1)) or len(range) == 1:
        #in both these cases, only one run is required, so pp not required
        if type(selection) == type({}) : selection = [selection]
        # put range into the correct format:
        if type(range) == type([]) or type(range) == type(arange(0,1)):
            range = range[0]
        for i in xrange(len(selection)):
            this_sel = deepcopy(selection[i])
            this_sel[variable] = range
            selections.append(this_sel)
        classes = make_detections(parameter,colour,cherrypick,index,sourcedata,selections,perfect_clusters=perfect_clusters,set_density=set_density)
        if len(classes) > 1: return classes,server
        else: return classes[0],server

# this is the case where (eg) one selection in, many out - defined by array of variable values...
# this routine does *not* handle many in, many out - ie len(selection) and len(range) > 1
    for entry in range:
        this_sel = deepcopy(selection)
        try:
            this_sel[variable] = entry
            this_sel = [this_sel]
        except KeyError:
            print '\t#Variable %s not in this selection!'
            return None,server
        selections.append(this_sel)
    args = [{'parameter':parameter},{'colour':colour},{'cherrypick':cherrypick},{'index':index},{'sourcedata':sourcedata},{'source[i]':None}]
    if pp: 
        if not perfect_clusters: args.append({'perfect_clusters':perfect_clusters})
        import signal
        classes,server = pp_run(selections,'make_detections',args,dependent_functions,modules,ncpus=ncpus,idle_delay=1200.0,server=server) # Behzad changed idle_delat from 30.0 to 1200.0
    else:
        for i in xrange(len(selections)):
            classes.append(make_detections(parameter,colour,cherrypick,index,sourcedata,selections[i],perfect_clusters=perfect_clusters))

    if len(classes) == 1: classes = classes[0]
    else:
        result = []
        for c in classes:
            result.append(c[0])
        classes = result
    for c in classes:
        c.clusters = []

    return classes,server

def make_mask(data,sel,plot=False,colour='gr'):
    passback = False
    mask = clusters.qselection_function('x',data,sel,plot=plot)
    return mask

def selection_fn(array):
    from tools import element_area
    observed = False
    grslope,grintercept,grwidth = array[0][0],array[0][1],array[0][2]
    rislope,riintercept,riwidth = array[1][0],array[1][1],array[1][2]
    izslope,izintercept,izwidth = array[2][0],array[2][1],array[2][2]
    zyslope,zyintercept,zywidth = array[3][0],array[3][1],array[3][2]
    if array[4][0] and array[4][1]:
        zmin,zmax,zeq = array[4][0],array[4][1],array[4][2]
    else:
        observed = True
        zmin,zmax,zeq = None,None,None

    
    xmin,xmax,ymin,ymax = array[5][0],array[5][1],array[5][2],array[5][3]
    rlimit = array[6]
    probthresh = array[7]
    smallest_cluster = array[8]
    rho_crit = array[9]
    sourcefile = array[10]
    area = array[11]
    halomass = array[12]
    selection = {'xmin':xmin,'xmax':xmax,'ymin':ymin,'ymax':ymax, \
                 'gr_slope':grslope, 'gr_intercept':grintercept,'gr_width':grwidth,\
                 'ri_slope':rislope, 'ri_intercept':riintercept,'ri_width':riwidth, \
                 'iz_intercept':izintercept,'iz_width':izwidth,'iz_slope':izslope,\
                 'zy_intercept':zyintercept,'zy_width':zywidth,'zy_slope':zyslope,\
                 'zeq':zeq,'rlimit':rlimit,'zmin':zmin,'zmax':zmax,'probthresh':probthresh,\
                 'smallest_cluster':smallest_cluster,'rho_crit':rho_crit,'sourcefile':sourcefile,\
                 'halomass':halomass}
    if observed:
        selection['observed_data'] = True
#    selection['detection_area'] = (degrees(selection['xmax']-selection['xmin'])*degrees(selection['ymax']-selection['ymin']))
    selection['detection_area'] = element_area(degrees(selection['xmin']),degrees(selection['xmax']),degrees(selection['ymin']),degrees(selection['ymax']))

    
    
#    selection['probthresh_function'] = 
    

    return selection


def cutData(data,colours,parameters,rlimit,ncpus=8):
    from numpy import arange,take,where,array
    import detections as classes
    from operator import countOf,indexOf,contains
    from copy import copy,deepcopy
    import gc
#    import fits_table as fits
    from lib.pp import pp as ppserver
    import time
    
    if ncpus > 0:
        pp = True
        ppservers = ()
        dependent_functions = (selection_fn,slice_data,make_mask,take,make_detections,classes)
        modules = ('clusters','numpy','parallel_qhull','operator','fits_table','operator','detections')
        ppthreshold = int((ncpus/4))
    else:
        pp = False
        ppthreshold = 1E10

    colour = False
    parameter = False
    cherrypick = True

    if type(data) != type({}):
        from dataIO import getData
        if data[1] == 'parameter': cherrypick = False
        data = getData(data[0])

    # determine the iteration - in colour space or parameters:
    col_query = [sum(colours[i][2]) > 0. for i in xrange(len(colours))]
    if countOf(col_query,True) == 1: colour = True
    if countOf(col_query,True) > 1:
        print '\tMore than one variation in colour space!'
        return None

    par_query = [parameters[2][i] != 0. for i in xrange(len(parameters[2]))]
    if countOf(par_query,True) == 1: parameter = True
    if countOf(par_query,True) > 1 and pp:
#        from detect import update
        catalogue = []
        catalogue.append(data['sourcefile'])
        catalogue.append('parameter')
        
        Names = ['probthresh','smallest_cluster','xc','yc','win','redshift','dredshift','rho_crit']
        index = where(par_query[:])[0][0]
        index2 = where(par_query[:])[0][1]
        print '\t#Exploring %s and %s in parameter space'%(Names[index],Names[index2])
        v1,v2,dv = parameters[0][index],parameters[1][index],parameters[2][index]
        w1,w2,dw = parameters[0][index2],parameters[1][index2],parameters[2][index2]
        v = arange(v1,v2+0.05*dv,dv)
        w = arange(w1,w2+0.05*dw,dw)
        print '\t#changing %s from %f to %f in increments %f:'%(Names[index],v1,v2,dv)
        print v
        print '\t#changing %s from %f to %f in increments %f:'%(Names[index2],w1,w2,dw)
        print w
        vars = [v,w]
        var_index = [index,index2]
        ############################## PARALLEL PYTHON INIT ################################
        ppservers = ()
        # have an array m*n. To reduce disk I/O and perfect cluster/class production, submit
        # min(m,n) jobs - ie fewest number of parallel submissions possible.
        optimal_var = vars[indexOf([len(l) for l in vars],min([len(l) for l in vars]))]
        optimal_index = var_index[indexOf([len(l) for l in vars],min([len(l) for l in vars]))]
        
        source = []
        # produce an array of the parameters...
        p = parameters
        for i in xrange(len(optimal_var)):
            start = []
            finish = []
            change = []
            #reproduce the range of parameters:
            for j in xrange(len(parameters[2])):
                start.append(p[0][j])
                finish.append(p[1][j])
                change.append(p[2][j])
            #now replace the key indices - ie, change the variable in var to change=0,
            #and the ith value from the vars array:
            start[optimal_index] = optimal_var[i]
            #finish tuple irrelevant, as code looks for change and start, but.....
            finish[optimal_index] = optimal_var[i]
            #ensure the change value is zero, so it doesn't iterate through this parameter:
            change[optimal_index] = 0.0
            #now tuple-ise:
            start = tuple(start)
            finish = tuple(finish)
            change = tuple(change)
            env = (start,finish,change)
            #add to the source:
            source.append(env)

        completed = []
        if pp: print '\t#Submitting %d colour slice productions to parallel job server'%(len(source))       
        args = [{'catalogue':catalogue},{'colours':colours},{'source[i]':None},{'rlimit':rlimit}]
        results = pp_run(source,'cutData',args,dependent_functions,modules,ncpus=ncpus)
        output = []
        for batch in results:
            for job in batch:
                output.append(job)
        results = output
        return results

    if colour and parameter:
        print '\t!!Variation in colour AND parameter space!'
        return []

    if not colour and not parameter:
        print '\t#No variations detected, producing one set of selected galaxies'
        index = -1
        selections = []
        for k in range(len(colours)):
                a = [colours[j][0] for j in range(len(colours))]
                for m in range(len(a)):
                    a[m] = [a[m][n] for n in range(len(a[m]))]

        z = (parameters[0][5]-parameters[0][6],parameters[0][5]+parameters[0][6],None)
        xmin,xmax = parameters[0][2]-parameters[0][4],parameters[0][2]+parameters[0][4]
        ymin,ymax = parameters[0][3]-parameters[0][4],parameters[0][3]+parameters[0][4]
        xy = [xmin,xmax,ymin,ymax]
        probthresh = parameters[0][0]
        smallest_cluster = parameters[0][1]
        rho_crit = parameters[0][7]
        a.append(z)
        a.append(xy)
        a.append(rlimit)
        a.append(probthresh)
        a.append(smallest_cluster)
        a.append(rho_crit)
        a.append(data['sourcefile'])
        a.append((degrees(parameters[0][4])*2)**2)
        a.append(parameters[0][8])
        sel = selection_fn(a)
        selections.append(sel)
                    
    if colour:
        cherrypick = True
        z = (parameters[0][5]-parameters[0][6],parameters[0][5]+parameters[0][6],None)
        xmin,xmax = parameters[0][2]-parameters[0][4],parameters[0][2]+parameters[0][4]
        ymin,ymax = parameters[0][3]-parameters[0][4],parameters[0][3]+parameters[0][4]
        xy = [xmin,xmax,ymin,ymax]
        probthresh = parameters[0][0]
        smallest_cluster = parameters[0][1]
        rho_crit = parameters[0][7]
        Names = ['gr','ri','iz','zy']
        Variables = ['slope','intercept','width']
        index = where(col_query[:])[0][0]
        index2 = where([colours[index][2][i] > 0 for i in xrange(len(colours[index][2]))])[0][0]
        print '\t#Exploring %s in %s colour space'%(Variables[index2],Names[index])
        print '\t#changing %s %s from %f to %f in increments %f'%(Names[index],Variables[index2],colours[index][0][index2],colours[index][1][index2],colours[index][2][index2])
        v1,v2,dv = colours[index][0][index2],colours[index][1][index2],colours[index][2][index2]
        array = arange(v1,v2+0.05*dv,dv)
        print array
        if len(array) > ppthreshold:
#            from detect import update
            catalogue = []
            catalogue.append(data['sourcefile'])
            catalogue.append('colour')
        ############################## PARALLEL PYTHON INIT ################################
            
            # min(m,n) jobs - ie fewest number of parallel submissions possible.
#            optimal_var = vars[indexOf([len(l) for l in vars],min([len(l) for l in vars]))]
#            optimal_index = var_index[indexOf([len(l) for l in vars],min([len(l) for l in vars]))]
        
            source = []
            # produce an array of the parameters...

            #which is the colour array that we are varying over?

            for i in xrange(len(array)):
                thisrun = list(deepcopy(colours))
                thisrun[index] = list(thisrun[index])
                thisrun[index][0] = list(thisrun[index][0])
                thisrun[index][1] = list(thisrun[index][1])
                thisrun[index][2] = list(thisrun[index][2])
                thisrun[index][0][index2] = array[i]
                thisrun[index][1][index2] = array[i]
                thisrun[index][2][index2] = 0.0
                source.append(thisrun)

            completed = []
            if pp: print '\t#Submitting %d colour slice productions to parallel job server'%(len(source))       
            
            args = [{'catalogue':catalogue},{'source[i]':None},{'parameter':parameters},{'rlimit':rlimit}]
            results = pp_run(source,'cutData',args,dependent_functions,modules,ncpus=ncpus)

            return results

    if parameter:
        cherrypick = False
        Names = ['probthresh','smallest_cluster','xc','yc','win','redshift','dredshift','rho_crit']
        index = where(par_query[:])[0][0]
#        index2 = where([colours[index][2][i] > 0 for i in range(len(colours[index][2]))])[0][0]
        print '\t#Exploring %s in parameter space'%(Names[index])
        v1,v2,dv = parameters[0][index],parameters[1][index],parameters[2][index]
        print '\t#changing %s from %f to %f in increments %f:'%(Names[index],v1,v2,dv)
        array = arange(v1,v2+0.05*dv,dv)
        print array


    if colour:
        selections = []
        for i in xrange(len(array)):
            for k in xrange(len(colours)):
                a = [colours[j][0] for j in xrange(len(colours))]
                for m in xrange(len(a)):
                    a[m] = [a[m][n] for n in xrange(len(a[m]))]
#            print '\t!!changing %s %s (value at a[%d][%d]=%f) to %f'%(Names[index],Variables[index2],index,index2,a[index][index2],array[i])
            a[index][index2] = array[i]
#            for m in range(len(a)):
#                print Names[m],a[m]
            a.append(z)
            a.append(xy)
            a.append(rlimit)
            a.append(probthresh)
            a.append(smallest_cluster)
            a.append(rho_crit)
            a.append(data['sourcefile'])
            a.append((degrees(parameters[0][4])*2)**2)
            a.append(parameters[0][8])
            sel = selection_fn(a)
            selections.append(sel)
    if parameter:
        selections = []
        for i in xrange(len(array)):
            for k in xrange(len(colours)):
                a = [colours[j][1] for j in xrange(len(colours))]
                for m in range(len(a)):
                    a[m] = [a[m][n] for n in xrange(len(a[m]))]


            z = (parameters[0][5]-parameters[0][6],parameters[0][5]+parameters[0][6],None)
            xmin,xmax = parameters[0][2]-parameters[0][4],parameters[0][2]+parameters[0][4]
            ymin,ymax = parameters[0][3]-parameters[0][4],parameters[0][3]+parameters[0][4]
            xy = [xmin,xmax,ymin,ymax]
            probthresh = parameters[0][0]
            smallest_cluster = parameters[0][1]
            rho_crit = parameters[0][7]
            if index == 0: probthresh = array[i]
            if index == 1: smallest_cluster = array[i]
            if index == 2: xmin,xmax = array[i]-parameters[0][4],array[i]+parameters[0][4]
            if index == 3: ymin,ymax = array[i]-parameters[0][4],array[i]+parameters[0][4]
            if index == 4:
                xmin,xmax = parameters[0][2]-array[i],parameters[0][2]+array[i]
                ymin,ymax = parameters[0][3]-array[i],parameters[0][3]+array[i]
                xy = [xmin,xmax,ymin,ymax]
            if index == 5: z = (array[i]-parameters[0][6],array[i]+parameters[0][6],None)
            if index == 6: z = (parameters[0][5]-array[i],parameters[0][5]+array[i],None)
            if index == 7: rho_crit = array[i]
            a.append(z)
            a.append(xy)
            a.append(rlimit)
            a.append(probthresh)
            a.append(smallest_cluster)
            a.append(rho_crit)
            a.append(data['sourcefile'])
            a.append((degrees(parameters[0][4])*2)**2)
            a.append(parameters[0][8])
            sel = selection_fn(a)
            selections.append(sel)

    print '\t#Produced selection functions'
    results = make_detections(parameter,colour,cherrypick,index,data,selections)
#    print '\tNew method complete'

    gc.collect()
    return results

def question(string,answer_type,escape=None,mixed_mode=True):
    """\
    question(string,answer_type,escape=None,mixed_mode=True) -> answer

    \t Asks the user a (blocking) question, requires answer to be of
    a specific type, or optionally an escape character.

        string     :  the question to be asked
        answer_type:  what format should the response be?   [string,int,float]
        escape     :  set an escape sequence that returns None  (e.g. vi-stype !Q)
        mixed_mode :  permit numbers in string       (e.g 'z=0.12')

    """

    import re
    if 1:
        def typeError(yes_type,no_type,escape=None):
            if escape != None:
                print """\tExpecting data type %s as input, got %s. Remember you can exit with '%s'"""%(answer_type,_type,escape)
            else:
                print '\tExpecting data type %s as input, got %s.'%(answer_type,_type)
            return


    import numpy
    types = {}
    types['string'] = [str,numpy.str]
    types['int'] = [numpy.int,numpy.int0,numpy.int8,numpy.int16,numpy.int32,numpy.int64,int]
    types['float'] = [numpy.float,numpy.float32,numpy.float64,numpy.float128,float]
    answer = None
    if answer_type not in types.keys():
        print 'Answer_type must be either: %s'%(types.keys())
    if string[-1] != ' ':
        string = string + ' > '
        
    while answer == None:

        try:
            answer = raw_input(string)
        except EOFError:
#            print 'Break detected'
            raise SystemExit('Break detected')

        if answer == escape:
            raise SystemExit('Break detected')


        if answer_type == 'int':
            if '.' in answer:
                _type = 'float'
                
                typeError(answer,_type,escape=escape); answer = None
            elif re.match("^[0-9]*$", answer):
                return int(answer)
        elif answer_type == 'float':
            if '.' in answer:
                if re.match("^[0-9]*$", ''.join(answer.split('.'))):
                    return float(answer)
                else:
                    _type = 'string'
                    typeError(answer,_type,escape=escape); answer = None
                    
        elif answer_type == 'string':
            nums = False
            letters = False
            for a in answer:
                if re.match("^[0-9]*$",a):
                    nums = True
                if re.match("^[A-Za-z]]*$",a):
                    letters = True

            if mixed_mode:
#                if (not letters) and nums:
#                    print 'No letters found in this string, only numbers! (but allow)?'
                return answer
                
            else:
                if letters and nums:
                    _type = 'string/numbers'
                    typeError(answer,_type,escape=escape); answer = None
                elif letters:
                    return answer
        

def setcfgProp(dataset,property,value,cfgloc='auto',category=None,addslash=False,verbose=True):
    lw = 21
    reserved = ['filename','gridfile','ra','dec']
    if cfgloc == 'auto':
        cfgloc = getcfgProp('orca','cfgpath',cfgloc='./')
    import os
    if cfgloc[-1] != '/':
        cfgloc = cfgloc+'/'
    filename = cfgloc+dataset+'.cfg'
    if not os.path.isfile(filename):
        if verbose: print 'File %s does not exist'%(filename)
        return None
    file = open(filename,'r')
    instream = file.readlines()
    file.close()
    outstream = []
    lines = iter(instream)
    replace = False
    
    for l in instream:
        key = l.split(':')[0]
        
        key0 = key
        if key0[0] == ' ':
            key0 = key0.split()[-1]
        elif key0[-1] == ' ':
            key0 = key0.split()[0]

#        if key[0] == ' ':
#            key = key.split()[-1]
        if key0 == property:
            replace = True
#            print repr(key)
    _thiscat = None
#    print replace
#    raise SystemExit(0)
    if replace:
        file = open(filename,'w')
        for l in instream:
            key = l.split(':')[0]
            key0 = key
            if key0[0] == ' ':
                key0 = key0.split()[-1]
            elif key0[-1] == ' ':
                key0 = key0.split()[0]
            
            if key0 == property:
                _l = '%s%s: %s\n'%(key,(lw-(len(key)))*' ',str(value))
                outstream.append(_l)
            else:
                outstream.append(l)
    else:
        file = open(filename,'w')
        added = False
        while 1:
            try:
                l = lines.next()

                if l == '':
                    outstream.append(l)
                elif l == '\n':
                    if (category == _thiscat) and (added == False):
                        _property = '     '+property
                        _l = '%s%s: %s\n'%(_property,(lw-(len(_property)))*' ',str(value))
                        outstream.append(_l)
                        added = True
                    else:
                        outstream.append(l)
                elif l[0] == '#':
                    outstream.append(l)
                elif l[0] != ' ':
                    #either a new category (no val, or #) or a generic value
                    
                    key = l.split(':')[0]
                    key0 = key.split(' ')[0]
                    val = l.split(':')[1]
                    if val[0] == ' ':
                        if '\n' in val:
                            val0 = val.split('\n')[0]
                        else:
                            val0 = val.split()
                    if 0:
#                    except:
                        print l
                        print key
                        print key0
#                        val = l.split(':')[1]
                        raise SystemExit('Error')

#                    if 'extras' in l:
#                        print repr(l)
#                        print repr(key)
#                        print repr(val)
#                        print repr(val0)
#                        print val0[0] == '#'
#                        raise SystemExit('Halt')
                    if (val == '\n') or (val0[0] == '#') or '#' in val:
#                        print 'possible category found: '+key
                        _thiscat = key
#                        print repr(key)
                    if val != '\n':
                        ## generic value! add data, append and then add this current line
                        if (added == False) and ((key0 in reserved) == False):
                            if category == None:
                                _l = '%s%s: %s\n'%(property,(lw-(len(property)))*' ',str(value))
                                outstream.append(_l)
                                added = True
                    outstream.append(l)
                elif l[0] == ' ':
                    #value under a category...
#                    print 'category value: '+l
                    if (category == _thiscat) and (added == False):
                        _property = property
                        npad = 0
                        while l[npad] == ' ':
                            npad += 1
                        _property = '%s%s'%(npad*' ',property)
                        _l = '%s%s: %s\n'%(_property,(lw-(len(_property)))*' ',str(value))
                        outstream.append(_l)
                        added = True

                    outstream.append(l)
                else:
                    outstream.append(l)
            except StopIteration:
                if added == False:
                    a = 1/0.
                    raise SystemExit('Key not added!!')
                break
        if added == True:
            if verbose:
                cat_txt = ''
                if category:
                    cat_txt = """ (under category '%s')"""%(category)
                if replace:
                    print """Value of key '%s' replaced by value: '%s'%s"""%(property,str(value),cat_txt)
                else:
                    print """Value of key '%s' replaced by value: '%s'%s"""%(property,str(value),cat_txt)
    file.writelines(outstream)
    file.close()
        
                

def getcfgProp(dataset,property,cfgloc='auto',addslash=False,verbose=True):
    if cfgloc == 'auto':
        cfgloc = getcfgProp('orca','cfgpath',cfgloc='./')
    import os
    if '#' in property:
        if verbose: print 'Comment line asked for!'
        return None
#    print cfgloc
#    print '<<%s>>'%(cfgloc[-1])
    if cfgloc[-1] != '/':
        cfgloc = cfgloc+'/'
    filename = cfgloc+dataset+'.cfg'
#    print filename
    if not os.path.isfile(filename):
        if verbose: print 'File %s does not exist'%(filename)
        return None
    file = open(filename)
    instream = file.readlines()
    result = None
    for line in instream:
#        print line
        if '#' in line and '#' == line[0]:
            continue
        identifier = line.split(':')[0]
        identifier = identifier.split('\t')[0]
        identifier = identifier.split(' ')
#        try:
#            identifier = identifier.split('\t')[0]
#        except:
#            pass
#        print identifier
        if property in identifier:
#            print line
            result = []
            
            line = line.split(':')[-1]
            line = line.split(' ')

            for point in line:
                if '\n' in point:
                    point = point.split('\n')[0]
                if point != '' and point != '\n':
                    result.append(point)
            if len(result) == 1:
                result = result[0]

            if addslash and result[-1] != '/':
                result += '/'


            return result

    if result == None:
        if verbose: print 'No match found'
    elif len(result) == 0:
        if verbose: print 'Match found, but no appropriate data extracted'


    return result

        

def make_cfg(no_clobber=False,indent=5):
    from tools import rows_cols
    import os
    global outstrings,entries,records
    
    trues = ['T','t','Y','y','True','true','Yes','yes','YES','TRUE']
    falses = ['F','f','N','n','False','false','No','no','NO','FALSE']


    def prepare(string,margin_len=21,term=':',indent=0):

        if indent != 0:
            _label = ''
            while(len(_label)) != indent:
                _label = ' '+_label
            string = _label+string

        if len(string)> margin_len:
            return string+'%s '%(term)

        while(len(string)) != margin_len:
            string = string + ' '
        string = string+'%s '%(term)
        return string

    def parse(string):
        global outstrings
        if string[-1:] != '\n':
            string = string +' \n'
        outstrings.append(string)
            
    def process_entry(string,label,noref=False,verbose=True,indent=0):
        global records,outstrings
        _string = ''
        data = []
        if string and string != '':
            for _s in string.split(' '):
                if _s == '':
                    continue
                try:
                    j = int(_s)
                    if j in records.keys() and noref==False:
                        _string = _string + '%s '%(records[j])
                        data.append(records[j])
                    else:
                        _string = _string + '%s '%(_s)
                        data.append(_s)
                except ValueError:
                    _string = _string + '%s '%(_s)
                    data.append(_s)
            

        if indent != 0:
            _label = ''
            while(len(_label)) != indent:
                _label = ' '+_label
            label = _label+label
        line = prepare(label)+_string;parse(line)
        if verbose:
            print line
        return data

    def fits_entries(indent=5):
        global entries
        
        print '\n'
        sourcecat = ["""Source catalogue columns:"""]
        c1 = '# Ref  '
        _c1 = prepare(c1,margin_len=len(c1),term='#')
        c2 = ' Entry '
        lmax = max(max([len(_s) for _s in entries])+3,len(c2))
        c2 = prepare(c2,margin_len=lmax,term='#')
        sourcecat.append(_c1+c2)

        _c1 = prepare('# ',margin_len=len(c1),term='#')
        _c2 = prepare(' ',margin_len=lmax,term='#')

        l2 = '#'
        while len(l2) < len(sourcecat[1])-1:
            l2 = l2 + '#'
        sourcecat.append(l2)
        sourcecat = ["""Source catalogue columns:""",l2]+sourcecat[1:]

        for i in xrange(len(entries)):
#            records[i+1] = tbl.names[i]

            _c1 = prepare('# %d'%(i+1),margin_len=len(c1),term='#')
            _c2 = prepare(' %s'%(tbl.names[i]),margin_len=lmax,term='#')

            sourcecat.append(_c1+_c2)
        sourcecat.append(l2)
        sourcecat.append(' ')

        for s in sourcecat:
            print s



    import os
    import numpy
    margin_len = 21        # ie ' '*21, then ':' 
    outstrings = []
    print 'ORCA cfg generator'
    cfgloc = getcfgProp('orca','cfgpath',cfgloc='./')
    print 'Will write .cfg file to %s   [to change this, quit and edit your orca.cfg file]'%(cfgloc)
    dataset = None
    while dataset in [None,'']:
        dataset = question('Please enter the name of this dataset [will store as <%s/dataset.cfg>]'%(cfgloc),'string',escape='!Q')


#        dataset = 'dataset';print 'REMOVE THIS AFTER TESTING'

#    dataset = '%s/%s.cfg'%(cfgloc,dataset)
    
    if dataset == None:
        return
    if cfgloc[-1] != '/':
        cfgloc = cfgloc+'/'
    if os.path.isfile(cfgloc+dataset+'.cfg') and no_clobber:
        print 'Sorry, cfg file for that dataset already exists. Either delete it, or re-run with no_clobber=False'
        return None
    cfgfile = cfgloc+dataset+'.cfg'

    sourcefile = 'None'
    while sourcefile == 'None' and not os.path.isfile(sourcefile):
        sourcefile = question('Please enter the absolute path of the source .fit[s] file','string',escape='!Q')
#        sourcefile = '/data/rw1/dmurphy/data/md01/sextractor_cat_merge_proc.fits';print 'REMOVE THIS AFTER TESTING'
#        sourcefile = '/export/disk/dmurphy/data/md01/sextractor_cat_merge_proc.fits';print 'REMOVE THIS AFTER TESTING'

        if sourcefile == None:
            return

        elif not os.path.isfile(sourcefile):
            print 'File %s does not exist, please re-enter'%(sourcefile)
            sourcefile = 'None'
    
    line = prepare('filename')+sourcefile
    parse(line)

    gridfile = question('Enter the absolute path of the grid file if it exists \n[Enter skips, option to generate available at end]','string',escape='!Q')
    if gridfile in [None,'']:
        process_entry('','grid',indent=0)
        gridfile = None        
    else:
        process_entry(gridfile,'grid',indent=0)

#    import fits_table as fits
    try:
        tbl = fits.fits_table(sourcefile)
    except:
        print 'Cannot read .fits file - is it corrupted?'
        return
    records = {}
    inv_records = {}
    entries = deepcopy(tbl.names)

    for i in xrange(len(entries)):
        records[i+1] = entries[i]
    for key in records.keys():
        inv_records[records[key]] = key

    fits_entries()

    ##########
    ra = question("""Enter the Ref number for the RA followed by 'r' if in radians (e.g 1 r):""",'string',escape='!Q')
    _ra = str(records[int(ra.split(' r')[0])])
    if 'r' in ra:
        _ra = _ra+' radians'
#    line = prepare('ra')+_ra;parse(line)


    ##########
    dec = question("""Enter the Ref number for the declination followed by 'r' if in radians (e.g 2 r):""",'string',escape='!Q')
    if 'r' in dec:
        dec = dec.split(' r')[0]+' radians'





    parse('positions:')
    ra = process_entry(ra,'ra',indent=5)[0]
    dec = process_entry(dec,'dec',indent=5)[0]
    
    ##########
    specz = question("""Enter the Ref number[s] for speczs [Enter skips]:""",'string',escape='!Q')
    specz = process_entry(specz,'spectroscopic_z',indent=indent)
    if len(specz) > 1:
        print specz
        specwt = question('Enter the statistical weights of these redshifts if required (e.g. 1 4 2) [Enter skips]','string',escape='!Q')
        if specwt:
            process_entry(specwt,'spectroscopic_wt',indent=indent,noref=True)
        else:
            line = prepare('spectroscopic_wt',indent=indent);parse(line)
            


    ##########
    photoz = question("""Enter the Ref number[s] for photozs [Enter skips]:""",'string',escape='!Q')
    photoz = process_entry(photoz,'photometric_z',indent=indent)
    if len(specz) > 1:
        print photoz
        photowt = question('Enter the statistical weights of these redshifts if required (e.g. 1 4 2) [Enter skips]','string',escape='!Q')
        if photowt:
            process_entry(photowt,'photometric_wt',indent=indent,noref=True)
        else:
            line = prepare('photometric_wt',indent=indent);parse(line)

    parse('bands:')
    _band_names = []
    data_bandnames = []
    col_ids = None
    band_col_ids = []
    while len(band_col_ids) < 2:
        band_col_ids = question("""Enter the Ref numbers for the bands, in the order ORCA will analyse (ie blue=>red) :""",'string',escape='!Q')
        band_col_ids = band_col_ids.split(' ')
        _bci = []
        for id in band_col_ids:
            if id != '':
                _bci.append(id)
        band_col_ids = _bci
        if len(band_col_ids) < 2:
            print 'ORCA requires at least two bands to operate'
    ##########
    band_names = [records[int(j)] for j in band_col_ids]
    print band_names
    user_band_names = []

    while len(band_col_ids) != len(user_band_names):
        user_band_names = question("""Enter the preferred names for these filters in the order above [Enter skips]':""",'string',escape='!Q')
#        print band_names, type(band_names)
        if user_band_names == None or user_band_names == '':
#            band_names = _bands
            break
        
        user_band_names = user_band_names.split(' ')
        
        if len(user_band_names) != len(band_col_ids):
            print 'Need to enter %d names'%(len(band_col_ids))

    if user_band_names:
        process_entry(' '.join(user_band_names),'band_names',indent=indent)
#        user_band_names = user_band_names.split(' ')
    else:
        process_entry(' '.join(band_names),'band_names',indent=indent)
        user_band_names = band_names
#        band_names = band_names.split(' ')
#        band_names = [records[int(_k)] for _k in __bands.split()]

    band_names = process_entry(' '.join(band_col_ids),'column_id',indent=indent)        
    band_lookup = {}          ## value is the fits column name
    inv_band_lookup = {}
    for i in xrange(len(band_col_ids)):
        band_lookup[user_band_names[i]] = band_names[i]
        inv_band_lookup[band_names[i]] = user_band_names[i]


#        line = prepare('photometric_wt',indent=indent);parse(line)    
    ##########
    ext = question("""Enter the Ref numbers for extinctions in these filters [Enter skips]':""",'string',escape='!Q')
    if ext:
        process_entry(ext,'extinction',indent=indent)
    else:
        line = prepare('extinction',indent=indent);parse(line)

    lims = []
    while (len(lims) != len(band_names)):
        lims = question("""Enter the flux limits for these filters (e.g. -99 24.1 23.7) \n[-99 = no limit, ? for interactive checker]':""",'string',escape='!Q')

        if lims == '?':
            _me = []
            magerrs = question("""Enter Ref numbers for filter flux errors [Enter = abort]':""",'string',escape='!Q')
            if magerrs not in [None,'']:
                magerrs = magerrs.split()
                for e in magerrs:
                    if int(e) in records.keys():
                        _me.append(int(e))
            magerrs = _me
            if len(magerrs) == len(band_names):
                if 'asts3' in socket.gethostname() == False:
                    plt.ion()
                fig = plt.figure()
                nrows,ncols = rows_cols(len(band_names))
                plots = iter(arange(len(band_names))+1)
                for i in xrange(len(band_names)):
                    b = tbl.get_column(band_names[i])
                    be = tbl.get_column(magerrs[i])
                    filter = (be < 0.18)
                    #w=0.152. half_w = 0.076. err_lim = half_w/50% = half_w/0.68 = 0.117
                    b = b[filter]
                    be = be[filter]
                    plt.subplot(nrows,ncols,plots.next())
                    plt.hexbin(b,be,bins='log',cmap=plt.cm.binary)
                    plt.xlabel(band_names[i])
                    plt.ylabel(band_names[i]+' error')
                    xlim = plt.xlim()
                    ylim = plt.ylim()
                    plt.plot([xlim[0],xlim[1]],[0.117,0.117],'r--')

            if 'asts3' in socket.gethostname():
                print 'This platform does not support matplotlib interactive.\nPlease make a note of data and close window'
                plt.show()
            else:
                plt.draw()
                plt.ion()
                plt.draw_if_interactive()
                plt.draw()
            lims = []


        elif len(lims.split()) != len(band_names):
            print lims
            print band_names
            print 'You must specify flux limits for each filter (or use -99 for each if no constraint)'
        else:
            lims = lims.split()
    process_entry(' '.join(lims),'limits',indent=indent)
    lim_lookup = {}
    for i in xrange(len(band_names)):
        lim_lookup[band_names[i]] = lims[i]
    print lim_lookup
        


    colours = []
    for i in xrange(len(user_band_names)):
        try:
            c = '%s-%s'%(user_band_names[i],user_band_names[i+1])
            colours.append(c)
        except:
            break
    cstring = ''
    for c in colours:
        cstring = cstring + '%s '%(c)
    cstring = cstring[:-1]
    bnstring = 'Band names: '
    for b in user_band_names:
        bnstring += '%s '%(b)
    bnstring = bnstring[:-1]
    print bnstring
    inp_colours = question("""Enter colours to be constructed from these bands (e.g. u-g g-r) [Enter for %s]':"""%(cstring),'string',escape='!Q')
    if inp_colours == None or inp_colours == '':
        inp_colours = cstring
        true_inp = ''
        for _c in inp_colours.split():

            _c1 = _c.split('-')[0]
            _c2 = _c.split('-')[1]
            
            if _c1 in user_band_names:
                b1 = _c1
            elif _c1 in band_names.keys():
                b1 = inv_band_lookup[_c.split('-')[0]]

            if _c2 in user_band_names:
                b2 = _c2
            elif _c2 in band_names.keys():
                b2 = inv_band_lookup[_c.split('-')[1]]



#            b2 = inv_band_lookup[_c.split('-')[1]]

#            b1 = band_lookup[_c.split('-')[0]]
#            b2 = band_lookup[_c.split('-')[1]]


            true_inp += '%s-%s '%(b1,b2)
        true_inp = true_inp[:-1]
        inp_colours = true_inp
    else:
        for _c in inp_colours.split():
            b1 = _c.split('-')[0]
            b2 = _c.split('-')[1]
            if b1 in inv_band_lookup.keys():
                true_inp += '%s-%s '%(b1,b2)        ##have entered originals, e.g u
            elif b1 in band_lookup.keys():
                b1 = band_lookup[_c.split('-')[0]]
                b2 = band_lookup[_c.split('-')[1]]
                true_inp += '%s-%s '%(b1,b2)        ##have entered preferred, e.g umag

        true_inp = true_inp[:-1] 
        inp_colours = true_inp 
    inp_colours = process_entry(inp_colours,'colours',indent=0)

    current_ranges = {}
    current_ranges['g-r'] = '[0.47,1.9]'
    current_ranges['r-i'] = '[-0.2,1.35]'
    current_ranges['i-z'] = '[-0.2,0.8]'
    current_ranges['z-y'] = '[-0.4,1.0]'

    cranges = None
    inter_plot = False
    while cranges in ['',None,'?']:
        cranges = question("""Colour ranges to search (e.g. [0.47,1.9] [-0.2,1.35]) \n[? for interactive checker, Enter for presets]':""",'string',escape='!Q')
        if cranges == None or cranges == '':
            cranges = ''
            for c in inp_colours:
                if c in current_ranges.keys():
                    cranges = cranges + '%s '%(current_ranges[c])
                elif '-' in c:
                    _this_range = None
                    _this_int_plot = False
                    while _this_range in ['',None,'?']:

                        _this_range = question("""%s not in presets: please enter explicitly (e.g. [0.47,1.9]) [? for interactive checker]"""%(c),'string',escape='!Q')
                        if _this_range == '?' and _this_int_plot == True:
                            print 'Already plotted this once - did you close it?!'
                            this_range = None
                        elif _this_range == '?':
                            if 'asts3' in socket.gethostname() == False:
                                plt.ion()
                            fig = plt.figure()
                            b1 = tbl.get_column(band_lookup[c.split('-')[0]])
                            b2 = tbl.get_column(band_lookup[c.split('-')[1]])
                            col = b1-b2
                            maglim = abs(float(lim_lookup[band_lookup[c.split('-')[1]]]))
                            filter = (col < 2) & (col > -2) & (b2 < maglim)
                            b2 = b2[filter]
                            col = col[filter]
                #            plt.plot(arange(5))
                            plt.hexbin(b2,col,cmap=plt.cm.jet,bins='log')
                            plt.ylim(-2,2)
                            plt.xlim(min(b2),maglim)
                            plt.ylabel(c)
                            plt.xlabel(c[0])
                            if 'asts3' in socket.gethostname() == False:
                                plt.draw()
                                _this_int_plot = True
                            else:
                                print 'This platform does not support matplotlib interactive.\nPlease make a note of data and close window'
                                plt.show()
                                
                        if _this_range not in ['?',None,'']:
                            if '[' in _this_range and ']' in _this_range and ',' in _this_range:
                                cmin,cmax = _this_range.split('[')[1].split(']')[0].split(',')
                                try:
                                    _cmin = float(cmin)
                                    _cmax = float(cmax)
                                except ValueError:
                                    _this_range = None
                                
                            else:
                                _this_range = None

                            if _this_range == None:
                                print 'Incorrect format: please enclose two numbers, separated by a comma, with square brackets. E.g.: [0.4,0.5]'


                        

                    cranges = cranges + '%s '  %(_this_range)  

        elif cranges == '?' and inter_plot == True:
            print 'Already plotted this once - did you close it?!'
            cranges = None
        elif cranges == '?':
    #        print 'Will do interactive thing here'
            from tools import rows_cols
            nrows,ncols = rows_cols(len(inp_colours))
            plots = iter(arange(len(inp_colours))+1)

            if 'asts3' in socket.gethostname() == False:
                plt.ion()
            fig = plt.figure()

            for c in inp_colours:
                b1 = tbl.get_column(band_lookup[c.split('-')[0]])
                b2 = tbl.get_column(band_lookup[c.split('-')[1]])
                col = b1-b2
                maglim = abs(float(lim_lookup[band_lookup[c.split('-')[1]]]))
                filter = (col < 2) & (col > -2) & (b2 < maglim)
                b2 = b2[filter]
                col = col[filter]
                plt.subplot(nrows,ncols,plots.next())
    #            plt.plot(arange(5))
                plt.hexbin(b2,col,cmap=plt.cm.jet,bins='log')
                  

                plt.ylim(-2,2)
                plt.xlim(min(b2),maglim)
                if c in current_ranges.keys():
                    xlim = plt.xlim()
                    ylim = plt.ylim()
                    _cr = current_ranges[c].replace(']','   ').replace('[','   ')
                    cmin,cmax = _cr.split(',')
                    plt.plot([xlim[0],xlim[1]],[float(cmin),float(cmin)],'w--')
                    plt.plot([xlim[0],xlim[1]],[float(cmax),float(cmax)],'w--')


                plt.ylabel(c)
                plt.xlabel(c[0])
                if 'asts3' in socket.gethostname() == False:
                    plt.draw()
    #        plt.subplots_adjust(wspace=0.3,right=0.97,left=0.10,bottom=0.07,top=0.97,hspace=0.3)    
            if 'asts3' in socket.gethostname() == False:
                plt.draw()
                plt.ion()
                plt.draw_if_interactive()
                plt.draw()
            else:
                print 'This platform does not support matplotlib interactive.\nPlease make a note of data and close window'
                plt.show()


            cranges = None

        else:
            cranges_bak = deepcopy(cranges)
            cranges = cranges.replace(']','   ').replace('[','   ').split()
            cout = ''
            if len(cranges) != len(inp_colours):
                print 'Incorrect number of ranges supplied - have %d, need %d'%(len(cranges),len(inp_colours))
                print cranges
                cranges = None
            else:
                for i in xrange(len(cranges)):
                    cr = cranges[i]                
                    if ',' in cr:
                        cmin,cmax = cr.split(',')
                        try:
                            _cmin = float(cmin)
                            _cmax = float(cmax)
                        except ValueError:
                            cr = None
                    else:
                        cr = None

                    if cr == None:
                        print 'Incorrect format for %s range: enclose two comma-separated numbers with square brackets. E.g.: [0.4,0.5]'%(inp_colours[i])
                        cranges = None
                        break
                if cranges != None:
                    #entry appears to check out OK
                    cranges = cranges_bak

    print cranges
    if type(cranges) == type([]):
        cranges = ' '.join(cranges)
    cranges = process_entry(cranges,'range',indent=0)

    ########################

    slopes = None
    while slopes == None:
        slopes = question("""Enter the slopes for these CMRs (e.g. -0.048 -0.017 -0.06 -0.10) [Enter = no slope data]:""",'string',escape='!Q')
        if slopes == '' or slopes == None:
            break
        else:
            if len(slopes.split()) != len(inp_colours):
                slopes = None
    if slopes:
        slopes = process_entry(slopes,'slope',indent=0)

    widths = None
    while widths == None:
        widths = question("""Enter the widths for these CMRs (e.g. 0.16 0.12 0.05 0.12) [Enter = no width data]':""",'string',escape='!Q')
        if widths == '' or widths == None:
            break
        else:
            if len(widths.split()) != len(inp_colours):
                widths = None
    if widths:
        widths = process_entry(widths,'widths',indent=0)

    fits_entries()
    objID = None
    while objID == None:
        objID = question("""Enter the Ref number for the object ID""",'string',escape='!Q')
        try:
            id = records[int(objID)]
        except:
            print 'please enter a valid ref number'
            objID = None

    process_entry(objID,'uid',indent=0)

    ext_spiel = """extras:              # anything here will be added to data['extras'] and eventually galaxy.extras   Format = <name>  : <data> || <name> : <data1> <data2> <data3>..."""
    line = prepare(ext_spiel,indent=0);parse(line)

    extras = question("""Enter the Ref number[s] for inherited extras [Enter skips]:""",'string',escape='!Q')
    if extras not in [None,'']:
        extras = extras.split()
        for e in extras:
            if int(e) in records.keys():
                extras_name = question("""Do you wish to define a name for column entry "%s" [Enter skips]:"""%(records[int(e)]),'string',escape='!Q')
                ename = records[int(e)]
                if extras_name not in [None,'']:
                    ename = extras_name
                process_entry(ename,records[int(e)],indent=indent)
                

    misc_spiel = """misc:                # set flags here, eg permit_dropouts : True     or angular_mask : True"""
    line = prepare(misc_spiel,indent=0);parse(line)
    miscs = {}
    miscs['colour_shuffle'] = 'Do you want to shuffle the galaxy colours (false detection estimator)?'
    miscs['maxBCG'] = 'Do you want to emulate the maxBCG input selection?'
    miscs['permit_dropouts'] = 'Does your catalogue contain drop-out entries (galaxies with flux not in all bands)?'
    for m in miscs.keys():
        ans = None
        while ans == None:
            ans = question(miscs[m]+' [T|F]','string',escape='!Q')
            if not ans in trues and not ans in falses:
                ans = None
            else:
                _ans = 'True'
                if ans in falses:
                    _ans = 'False'
                process_entry(_ans,m,indent=indent)

#    print ''.join(outstrings)
    file = open(cfgfile,'w')
    file.writelines(outstrings)
    file.close()

    if gridfile != None:
        print 'Written to %s'%(cfgfile)
        return

    while gridfile == None:
        gridfile = question('Would you now like to create a grid file and estimate the boundary of the data? [y|n]','string',escape='!Q')
        if gridfile in trues:
            print 'Calling ORCA Define ARbitrary Boundary (./darb_data)'
            os.system('./darb_data')
            is_ok = None
            while is_ok == None:
                darb_out = sourcefile.split('.fit')[0]+'.grid'
                darb_args0 = '-dataset %s -store True -src_override %s -ra %s -dec %s -outname %s'%(dataset,sourcefile,ra,dec,darb_out)
                print 'The following will be added to the argument list, please do not add them below:\n** %s **'%(darb_args0)
                darb_args = question('Enter optional args based on above (e.g.  -sf 70 -plot True) [Enter skips]','string',escape='!Q')
                if darb_args in [None,'']:
                    break
                forbidden = ['dataset','store','src_override','outname','ra','dec']
                for fb in forbidden:
                    bail = False
                    if fb in darb_args:
                        print 'Please do not override the %s argument'%(fb)
                    if bail:
                        continue
                string = './darb_data %s %s'%(darb_args0,darb_args)
                print string
                os.system(string)
                is_ok = question('Did DARB run as expected? [y|n, Enter to abort]','string',escape='!Q')
                if is_ok in trues:
                    break
                elif is_ok in falses:
                    is_ok = None
                elif is_ok in [None,'']:
                    darb_out = ' \n'
                    break
                

            print 'Should be done here'
            for i in xrange(len(outstrings)):
                line = outstrings[i]
                if 'grid           ' in line:
                    outstrings[i] = line.replace(' \n',darb_out+' \n')
                    break
            file = open(cfgfile,'w')
            file.writelines(outstrings)
            file.close()



        elif gridfile in falses:
            gridfile = 'n'
            


    print 'Written to %s'%(cfgfile)

    return



def getData(catalogue,plot=False,colour='gr',interact=False,dataset='s82',args=None,obs=False,wrap=True,extinction_correction=True,selection=None,argv=None,subset_index=None,phot_err=None,maxBCG_mode=True,grid=None,tbl_data=None,make_subs_bailsafe=False,stripe=None,no_constraint=False,raise_on_fail=False,cfg_linkage=False,md5=False):
#    import fits_table as fits
    import numpy
    import detections as classes
    from tools import element_area

    export = '/export/disk/dmurphy/'

    update = False
    sub_check = True
    sub_write = True

    if selection != None and selection['detection'].has_key('stripe') and (stripe == None):
        stripe = selection['detection']['stripe']

    print '#getData():'
    if type(catalogue) != type(''):
        allcats = catalogue
        catalogue = allcats[dataset]

    if selection!= None:
        if selection['detection'].has_key('subset_index') and selection['detection']['subset_index'] != None:
            print '#getData(): Auto-discovered subset index %d'%(selection['detection']['subset_index'])
            subset_index = selection['detection']['subset_index']
        
    search_colours = ['gr','ri','iz','zy']

    catalogue = None
    cfg0 = dataset+'.cfg'

    _cfg = getcfgProp('orca','cfgpath',cfgloc='./')
    if _cfg[-1] != '/':
        _cfg += '/'
    cfg = _cfg+'%s'%(cfg0)
    cfg2 = dataset+'.cfg'
    cfg3 = export+'data/%s/'%(dataset)+cfg0
    cfg4 = '/data/rw1/dmurphy/r3/PanSTARRS/cfg/%s'%(cfg0)
    

    x,y,umag,gmag,rmag,imag,zmag,ymag,specz,photoz,id = None,None,None,None,None,None,None,None,None,None,None
    vx,vy,vz,stellar_mass,halo_mass,halo_id,halo_idrep,halo_jm,halo_jmtree,halo_zcos= None,None,None,None,None,None,None,None,None,None

    uid = None
    uid_name = ''

    if not os.path.isfile(cfg) and not os.path.isfile(cfg2) and not os.path.isfile(cfg3) and not os.path.isfile(cfg4):
        orcaloc = getcfgProp('orca','cfgpath',cfgloc='./')
        print '\t#########################################################################'
        print '\t\n\n\n\n\t\tNo cfg file found for dataset <<%s>>\n\t\t[orca.cfg pointing to %s]\n\n\n'%(dataset,orcaloc)
        print '\t#########################################################################'
        raise SystemExit(0)

                        
    if not os.path.isfile(cfg):
        cfg = cfg2
    if not os.path.isfile(cfg):
        cfg = cfg3
    if not os.path.isfile(cfg):
        cfg = cfg4
    if not os.path.isfile(cfg):
        raise getDataError('#getData(): Cannot find .cfg file - need input!')

    _cfg = open(cfg)
    print '#getData(): reading configuration data from %s'%(cfg)
    instream = _cfg.readlines()
    line0 = instream[0]

    line_split = line0.split()
    if '*' in line_split[-1]:
        line0 = ' '.join(line0.split()[:-1])
        sub_check = False
        sub_write = False
        
    catalogue = line0.split()[-1]
    catalogue2 = export+'data/%s/'%(dataset)+line0.split()[-1]

    if os.path.isfile(catalogue) == False and os.path.isfile(catalogue2) == False:
        raise getDataError('#getData(): Catalogue file(s) missing. ORCA looked for:\n%s\n%s'%(catalogue,catalogue2))

    if os.path.isfile(catalogue2):
        catalogue = catalogue2

    
    cat_extra = catalogue.split('.fits')[0]+'_extra'+'.fits'
    if os.path.isfile(cat_extra):
        catalogue = cat_extra

    if 'dr7' in catalogue and 'dr7full' not in catalogue:
        cat_extra = catalogue.split('.fits')[0]+'_maxbcg'+'.fits'
        if os.path.isfile(cat_extra):
            if selection != None:
                if selection.has_key('overrides'):
                    if selection['overrides'].has_key('no_maxBCG'):
                        maxBCG_mode = False
            elif maxBCG_mode == False:
                maxBCG_mode = False
            else:
                print '#getData(): Found maxBCG-selected catalogue, overriding'
                catalogue = cat_extra
                maxBCG_mode = False
        else:
            if selection != None:
                if selection.has_key('overrides'):
                    if selection['overrides'].has_key('no_maxBCG'):
                        maxBCG_mode = False
            

    else:
        maxBCG_mode = False



    data = {}
    protected = {}
    filter_widths = {}
    filter_slopes = {}
    extinctions = {}
    permit_dropouts = False
    angular_mask = False
    colour_shuffle = False
    loc_shuffle = False
    gridfile = False
    offset_threshold = None
    cfg_cuts = {}
                
    data['mock_data'] = {'halo_mass':halo_mass,'halo_id':halo_id,'halo_jm':halo_jm,'halo_jmtree':halo_jmtree,'halo_zcos':halo_zcos,'halo_idrep':halo_idrep}
    data['spatial'] = {'x':x,'y':y,'photoz':photoz,'specz':specz}
    data['magnitudes'] = {}
    data['colours'] = {}
    data['extras'] = {'id':id,'vx':vx,'vy':vy,'vz':vz,'stellar_mass':stellar_mass}

    extra_data_input = False
    override_data_input = False
    overrides = {}
    crange = {}



    if 'fit' in catalogue:
#        import fits_table as fits
        input_name = catalogue
        new_name = ''

        stripe_str = ''
        if stripe:
            stripe_str = str(stripe)
            
            alt_name = catalogue.split('.fit')[0]+'%s.fit'%(stripe_str)
            if not ('%s.fit'%(stripe_str) in catalogue) and os.path.isfile(alt_name):
                input_name = alt_name

        if subset_index != None and sub_check:
            alt_name = catalogue.split('.fit')[0]+'%s.fit'%(stripe_str)

            if not ('%s.fit'%(stripe_str) in catalogue):
                new_name = catalogue.split('.fit')[0]+'%s_%d.fit'%(stripe_str,subset_index)
            else:
                new_name = catalogue.split('.fit')[0]+'_%d.fit'%(subset_index)

            if os.path.isfile(new_name):
                if make_subs_bailsafe:
                    print '#getData(): Subset data detected for sub_id=%d, bailing'%(subset_index)
                    return None,None,tbl_data,None,None
                print '#getData(): Using subset file %s'%(new_name)
                input_name = new_name
        
        print '#getData(): Data source is: %s'%(input_name)
#        raise SystemExit(0)


        if tbl_data and (os.path.isfile(new_name) == False):
            print '#getData(): Using pre-read data file'
            tbl = tbl_data

        else:
            print '\t#Reading input data from %s'%(input_name)
            tbl = fits.fits_table(input_name)            


        if not os.path.isfile(cfg):
            cfg = cfg2
        if not os.path.isfile(cfg):
            cfg = cfg3
        if not os.path.isfile(cfg):
            raise getDataError('Cannot find .cfg file - need input!')

        _cfg = open(cfg)
        instream = _cfg.readlines()
        line_numbers = {}
        for line in instream:
            if '#' in line and '#' == line[0]:
                    continue

            if 'grid' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers['grid'] = indexOf(instream,line)
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        gridfile = point
                if (gridfile != False) and os.path.isfile(gridfile):
                    print '#getData(): Will use gridfile %s'%(gridfile)
                elif gridfile == False:
                    pass
                else:
                    raise getDataError('#getData(): Specified grid data file does not exist - need input! \n\n\t\t(%s)'%(gridfile))

############ POSITIONS ###

            elif 'ra ' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers['ra'] = indexOf(instream,line)
                ra = ''
                rad = False
                if 'radians' in line:
                    #means that the input is already in radians
                    rad = True
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        if rad:
                            point = point.split()[0]
                            try:
                                data['spatial']['x'] = tbl.get_column(point)
                                cfg_cuts['x_rad'] = point
                            except AttributeError:
                                print '#getData(): Cannot read data - this might be an empty dataset - returning'
                                core = classes.detections(None,None,dummy=True)
                                core.selection = {}
                                core.selection['detection'] = {'empty':True}
                                core.selection['limits'] = {}
                                core.clusters = []
                                return None,None,None,None,core

                            ra_name = point
                            break
                        else:
                            try:
                                #source data is in degrees
                                data['spatial']['x'] = numpy.radians(tbl.get_column(point))
                                cfg_cuts['x_deg'] = point
                            except AttributeError:
                                print '#getData(): Cannot read data - this might be an empty dataset - returning'
                                core = classes.detections(None,None,dummy=True)
                                core.selection = {}
                                core.selection['detection'] = {'empty':True}
                                core.selection['limits'] = {}
                                core.clusters = []
                                return None,None,None,None,core

                            ra_name = point
                            break


            elif 'dec' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers['dec'] = indexOf(instream,line)
                dec = ''
                rad = False
                if 'radians' in line:
                    #means that the input is already in radians
                    rad = True
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        if rad:
                            point = point.split()[0]
                            data['spatial']['y'] = tbl.get_column(point)
                            cfg_cuts['y_rad'] = point
                            dec_name = point
                            break
                        else:
                            data['spatial']['y'] = numpy.radians(tbl.get_column(point))
                            cfg_cuts['y_deg'] = point
                            dec_name = point
                            break

            elif 'spectroscopic_z' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers['spectroscopic_z'] = indexOf(instream,line)
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        data['spatial']['specz'] = tbl.get_column(point)
                        cfg_cuts['specz'] = point
                        sz_name = point

            elif 'photometric_z' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers['photometric_z'] = indexOf(instream,line)
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        data['spatial']['photoz'] = tbl.get_column(point)
                        cfg_cuts['photoz'] = point
                        pz_name = point
############ BANDS & COLOURS ###

            elif 'band_names' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers['band_names'] = indexOf(instream,line)
                band_names = []
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        band_names.append(point)

                if len(band_names) == 0:
                    raise getDataError('#getData(): No band_names in .cfg file - need more input!') 
                    
                for name in band_names:
                    data['magnitudes'][name] = []

            elif 'column_id' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers['column_id'] = indexOf(instream,line)
                cid = []
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        cid.append(point)


                retain = numpy.array(numpy.zeros(len(data['spatial']['x'])),dtype='bool')
                if len(cid) == 0:
                   raise getDataError('#getData(): No coloumn_id names in .cfg file - need more input!') 

                for i in xrange(len(band_names)):
                    data['magnitudes'][band_names[i]] = tbl.get_column(cid[i])
                    cfg_cuts[band_names[i]] = cid[i]

                    protected[band_names[i]] = numpy.array(numpy.zeros(len(data['spatial']['x'])),dtype='bool')
                
            elif 'extinction' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers['extinction'] = indexOf(instream,line)
                eid = []
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        eid.append(point)
                src_magnitudes = {}

                if len(eid) == len(band_names):
                    print 'Applying extinction corrections....'
                    for i in xrange(len(band_names)):
                        src_magnitudes[band_names[i]] = data['magnitudes'][band_names[i]]
                        extinction = numpy.array(tbl.get_column(eid[i]))
                        extinctions[eid[i]] = extinction
                        cfg_cuts[eid[i]] = eid[i]
                        data['magnitudes'][band_names[i]] = data['magnitudes'][band_names[i]] - extinction
                elif len(eid) > 0:
                    raise getDataError('#getData(): Have %d but need %d extinction tags in .cfg file - need more input!'%(len(eid),len(band_names))) 
                elif len(eid) == 0:
                    for i in xrange(len(band_names)):
                        src_magnitudes[band_names[i]] = data['magnitudes'][band_names[i]]

                    

            elif 'limits' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers['limits'] = indexOf(instream,line)
                limits = []
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        limits.append(float(point))
                
                lims = {}
                if len(limits) == len(band_names):
                    for i in xrange(len(band_names)):
                        if limits[i] != float(-99):
                            lims[band_names[i]] = limits[i]
                        elif limits[i] == float(-99):
                            lims[band_names[i]] = 40.0
                else:
                    raise getDataError('#getData(): Have %d but need %d magnitude limits in .cfg file - need more input!'%(len(limits),len(band_names))) 
                limits = lims

            elif 'colours' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers['colours'] = indexOf(instream,line)
                col_names = []
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
#                        point = ''.join(point.split('-'))
                        col_names.append(point)
                colours = {}
                new_col_names = []
                original_colours = col_names
                for i in xrange(len(col_names)):
                    col = col_names[i]
                    m1_name = col.split('-')[0]
                    m2_name = col.split('-')[-1]
                    m1 = data['magnitudes'][col.split('-')[0]]
                    m2 = data['magnitudes'][col.split('-')[-1]]
                    col = ''.join(col.split('-'))
                    data['colours'][col] = m1-m2
                    
                    drop_outs = (m1 < 0) | (m2 < 0) | (m1 > 40.0) | (m2 > 40.0)
                    data['colours'][col][drop_outs] = 680008

                    keepers = (drop_outs == False)
                    protected[m1_name][keepers] = 1
                    protected[m2_name][keepers] = 1
                    retain[keepers] = 1

                    new_col_names.append(col)
                col_names = new_col_names
                sequence = col_names
            
              

            elif 'range' in line and ((True in [override_data_input,extra_data_input]) == False):
                try:
                    dc = args['dc']
                except:
                    dc = 0.02
                line_numbers['range'] = indexOf(instream,line)
                ranges = []
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        ranges.append(point)
                if len(ranges) == len(col_names):
                    for i in xrange(len(col_names)):
                        col = ''.join(col_names[i].split('-'))
                        lims = ranges[i].split('[')[-1]
                        lims = lims.split(']')[0]
                        _min,_max = lims.split(',')
                        _min = float(_min)
                        _max = float(_max)+(dc/2.0)
                        crange[col] = [_min,_max,dc]
                
                elif len(ranges) > 0:
                    raise getDataError('#getData(): Have %d but need %d colour ranges in .cfg file - need more input!'%(len(ranges),len(col_names))) 


            elif 'slope' in line and ((True in [override_data_input,extra_data_input]) == False):
                _ltmp = line
                line_numbers['limits'] = indexOf(instream,line)
                slopes = []
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        slopes.append(float(point))
                
                slps = {}
                if len(slopes) == len(col_names):
                    for i in xrange(len(col_names)):
                        if slopes[i] != float(-99):
                            slps[col_names[i]] = float(slopes[i])
                        elif slopes[i] == float(-99):
                            slps[col_names[i]] = 0.0
                elif len(ranges) > 0:
                    raise getDataError('#getData(): Have %d but need %d slopes in .cfg file - need more input!'%(len(slopes),len(col_names))) 

                filter_slopes = slps

            elif 'width' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers['limits'] = indexOf(instream,line)
                widths = []
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        widths.append(float(point))
                
                wdths = {}
                if len(widths) == len(col_names):
                    for i in xrange(len(col_names)):
                        if widths[i] != float(-99):
                            wdths[col_names[i]] = float(widths[i])
                        elif widths[i] == float(-99):
                            wdths[col_names[i]] = 100.0
                elif len(ranges) > 0:
                    raise getDataError('#getData(): Have %d but need %d widths in .cfg file - need more input!'%(len(widths),len(col_names))) 

                filter_widths = wdths




            elif 'uid' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers['uid'] = indexOf(instream,line)
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        data['extras']['id'] = tbl.get_column(point)
                        uid_name = point
############ MOCK DATA ###

            elif 'halomass' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers['halomass'] = indexOf(instream,line)
                halo_mass = ''
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        data['mock_data']['halo_mass'] = tbl.get_column(point)
                        cfg_cuts['halo_mass'] = point
                        if min(data['mock_data']['halo_mass']) > 50:
                            data['mock_data']['halo_mass'] = numpy.log10(data['mock_data']['halo_mass'])
                            
            elif 'halo_id ' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers[''] = indexOf(instream,line)
                halo_id = ''
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        data['mock_data']['halo_id'] = tbl.get_column(point)
                        cfg_cuts['halo_id'] = point
            elif 'halo_jm ' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers[''] = indexOf(instream,line)
                halo_jm = ''
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        data['mock_data']['halo_jm'] = tbl.get_column(point)
                        cfg_cuts['halo_jm'] = point
            elif 'halo_jmtree' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers[''] = indexOf(instream,line)
                halo_jmtree = ''
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        data['mock_data']['halo_jmtree'] = tbl.get_column(point)
                        cfg_cuts['halo_jmtree'] = point
            elif 'halo_zcos' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers[''] = indexOf(instream,line)
                halo_zcos = ''
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        data['mock_data']['halo_zcos'] = tbl.get_column(point)
                        cfg_cuts['halo_zcos'] = point
            elif 'halo_idrep' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers[''] = indexOf(instream,line)
                halo_idrep = ''
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        data['mock_data']['halo_idrep'] = tbl.get_column(point)
                        cfg_cuts['halo_idrep'] = point

############ EXTRAS ###

            elif 'velocity' in line:
                line_numbers['velocity'] = indexOf(instream,line)
                velocities = []
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        velocities.append(point)
                if len(velocities) > 0:
                    vels = iter(velocities)
                    _v = vels.next();data['extras']['vx'] = tbl.get_column(_v)
                    cfg_cuts['vx'] = _v
                    _v = vels.next();data['extras']['vy'] = tbl.get_column(_v)
                    cfg_cuts['vy'] = _v
                    _v = vels.next();data['extras']['vz'] = tbl.get_column(_v)
                    cfg_cuts['vz'] = _v



            elif 'stellar_mass' in line:
                line_numbers['stellar_mass'] = indexOf(instream,line)
                stellar_mass = ''
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        data['extras']['stellar_mass'] = tbl.get_column(point)
                        cfg_cuts['stellar_mass'] = point

            elif 'overrides' in line:
                override_data_input = True
                extra_data_input = False

            elif 'extras' in line:
                extra_data_input = True

            elif 'misc' in line:
                extra_data_input = False
                override_data_input = False


            elif extra_data_input or override_data_input:
                #means we have something not explictly defined, but will add for inheritance
                if '#' in line:
                    continue
                line_in = line
                identifier = line.split(':')[0]
                identifier = identifier.split(' ')
                id = ''
                for i in identifier:
                    if i != '':
                        id = i
                        break
                identifier = i
                entry = line.split(':')[-1]
                line = entry.split(' ')
                entries = []
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if point != '' and point != '\n':
                        entries.append(point)
                if len(entries) == 1:
                    if extra_data_input:
                        try:
                            data['extras'][identifier] = tbl.get_column(entries[0])
                            cfg_cuts[identifier] = entries[0]
                        except KeyError:
                            if raise_on_fail:
                                raise getDataError('#getData(): Cannot find <<%s>> in source data file - need better input!'%(identifier))
                            else:
                                print '\tWARNING: Could not find <<%s>> in source data - please change entry in .cfg file (%s)'%(identifier,cfg)
                    elif override_data_input:
                        _val = entries[0]
                        if _val in ['True', 'False']:
                            if _val == 'True':
                                _val = True
                            else:
                                _val = False
                        else:
                            #is string?
                            try:
                                _val_flt = float(_val)
                                try:
                                    _val_int = int(_val)
                                    #final check..
                                    if float(_val_int) == _val_flt:
                                        _val = _val_int
                                except ValueError:
                                    #float then...
                                    _val = _val_flt
                            except ValueError:
                                #_val is likely a string - leave it as it is..
                                pass
                        
                        overrides[identifier] = _val
                    else:
                        raise SystemExit('Neither extra nor override input....?')
                        
                elif len(entries) > 1:
                    for e in entries:
                        if extra_data_input:
                            try:
                                data['extras'][e] = tbl.get_column(e)
                                cfg_cuts[e] = e
                            except KeyError:
                                if raise_on_fail:
                                    raise getDataError('Cannot find <<%s>> in source data file - need better input!'%(identifier))
                                else:
                                    pass
                        elif override_data_input:
                            overrides[e] = e
                                    
                
                line_numbers[identifier] = indexOf(instream,line_in)


############ MISC ###

            elif 'multi_survey' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers['multi_survey'] = indexOf(instream,line)
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if '\t' in point:
                        point = point.split('\t')[0]
                    if point != '' and point != '\n':
                        try:
                            if point == 'False':
                                permit_dropouts = False
                            elif point == 'True':
                                permit_dropouts = True
                        except:
                            pass
            elif 'permit_dropouts' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers['permit_dropouts'] = indexOf(instream,line)
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if '\t' in point:
                        point = point.split('\t')[0]
                    if point != '' and point != '\n':
                        try:
                            if point == 'False':
                                permit_dropouts = False
                            elif point == 'True':
                                permit_dropouts = True
                        except:
                            pass


                        
            elif 'maxBCG' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers['maxBCG'] = indexOf(instream,line)
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if '\t' in point:
                        point = point.split('\t')[0]
                    if point != '' and point != '\n':
                        try:
                            if point == 'False':
                                maxBCG_mode = False
                            elif point == 'True':
                                maxBCG_mode = True
                        except:
                            pass

            elif 'angular_mask' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers['angular_mask'] = indexOf(instream,line)
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if '\t' in point:
                        point = point.split('\t')[0]
                    if point != '' and point != '\n':
                        try:
                            if point == 'False':
                                angular_mask = False
                            elif point == 'True':
                                angular_mask = True
                        except:
                            pass

            elif 'colour_shuffle' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers['colour_shuffle'] = indexOf(instream,line)
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if '\t' in point:
                        point = point.split('\t')[0]
                    if point != '' and point != '\n':
                        try:
                            if point == 'False':
                                colour_shuffle = False
                            elif point == 'True' and '_shuffle' not in input_name:
                                colour_shuffle = True
                        except:
                            pass

            elif 'loc_shuffle' in line and ((True in [override_data_input,extra_data_input]) == False):
                line_numbers['loc_shuffle'] = indexOf(instream,line)
                line = line.split(':')[-1]
                line = line.split(' ')
                for point in line:
                    if '\n' in point:
                        point = point.split('\n')[0]
                    if '\t' in point:
                        point = point.split('\t')[0]
                    if point != '' and point != '\n':
                        try:
                            if point == 'False':
                                loc_shuffle = False
                            elif point == 'True' and '_shuffle' not in input_name:
                                loc_shuffle = True
                        except:
                            pass

############ END ###

        if cfg_linkage:
            return cfg_cuts
        

        if gridfile and (not grid):
            from dataIO import load
            grid = load(gridfile)
            if stripe:
                grid = grid[stripe]

        if grid:
            if selection != None and selection['detection'].has_key('stripe') and (stripe == None):
                stripe = selection['detection']['stripe']

            if stripe == 1.142:
                pass
                if stripe in grid.keys():
                    grid = grid[stripe]
                    user_grid = grid
                else:
                    keys = grid.keys()
                    sid = [grid[k]['sub_id'] for k in keys]
                    if subset_index in sid:
                        grid = grid[keys[indexOf(sid,subset_index)]]
                        user_grid = grid

        if wrap and (grid != None or (selection != None and selection.has_key('grid'))) == False:
            from math import pi
            xform = (data['spatial']['x'] > 2.0)
            data['spatial']['x'][xform] = data['spatial']['x'][xform] - 2.0*pi

        elif stripe != None and ((grid !=None)) or (selection and selection['detection'].has_key('stripe_wrap')):
            stripe_wrap = None
            if (grid !=None):
                if stripe < 1:
                    print '#getData(): Removing superflous outer shell in grid for a non-defined stripe'
                    grid = grid[stripe]


                keys = grid.keys()
                sid = [grid[k]['sub_id'] for k in keys]
                if subset_index in sid:
                    _grid = grid[keys[indexOf(sid,subset_index)]]
                if _grid.has_key('wrap'):
                    stripe_wrap = _grid['wrap']

            elif (selection['detection'].has_key('stripe_wrap')):
                stripe_wrap = selection['grid_wrap']

            if stripe_wrap:
                #### stripe wrap will **always** displace by 2pi now....
                from numpy import radians
                offset_threshold = radians(stripe_wrap[0])
                xform = (data['spatial']['x'] < radians(stripe_wrap[0]))
                test = (data['spatial']['x'] > radians(360.0))
                data['spatial']['x'][xform] = data['spatial']['x'][xform] + radians(stripe_wrap[1])
                if True in xform:
                    xform = numpy.ones(len(xform),dtype='int')*xform
                elif True in test:
                    xform = numpy.ones(len(xform),dtype='int')*test
                    print '#getData(): Data was read in pre-wrapped, adding offsets to sources > RA=360'
#                xform = numpy.ones(len(xform))[xform]
                data['extras']['offset'] = numpy.array(xform,dtype='bool')


        else:
            from math import pi
            
            string =  """Keep an eye on this - might need to do something for the straddle point, ie

            xform = (data['spatial']['x'] < 0.0)
            data['spatial']['x'][xform] = data['spatial']['x'][xform] + 2.0*pi

                         """

            if min(data['spatial']['x']) < 0.0 and max(data['spatial']['x']) < 0.0 and '82' in dataset:
                data['spatial']['x'] = data['spatial']['x'] + 2.0*pi
                
                
        
        constraint = numpy.array(numpy.ones(len(data['spatial']['x'])),dtype='bool')

        if colour_shuffle and os.path.isfile(input_name.split('.fits')[0]+'_shuffle'+'.fits') == False and '_shuffle' not in input_name:
            import numpy
            update = True
            new_name = input_name.split('.fits')[0]+'_shuffle'+'.fits'
            indices = arange(len(data['spatial']['x']))
            numpy.random.shuffle(indices)
            print '#getData(): Now shuffling colours'

            outputs = {}
            ext_outs = {}
            for _band in data['magnitudes'].keys():
                outputs[_band] = []
            for _key in extinctions.keys():
                ext_outs[_key] = []
                
            for i in xrange(len(data['spatial']['x'])):
                for _band in outputs.keys():
                    outputs[_band].append(data['magnitudes'][_band][indices[i]])
                for _key in ext_outs.keys():
                    ext_outs[_key].append(extinctions[_key][indices[i]])
                    
            for _band in outputs.keys():
                data['magnitudes'][_band] = numpy.array(outputs[_band])
                src_magnitudes[_band] = numpy.array(outputs[_band])


            for _key in extinctions.keys():
                extinctions[_key] = numpy.array(ext_outs[_key])

            print '#getData(): Colours shuffled'


        if loc_shuffle and os.path.isfile(input_name.split('.fits')[0]+'_shuffle'+'.fits') == False and '_shuffle' not in input_name:
            import numpy
            update = True
            new_name = input_name.split('.fits')[0]+'_shuffle'+'.fits'
            indices = arange(len(data['spatial']['x']))
            numpy.random.shuffle(indices)
            print '#getData(): Now shuffling galaxy locations'

            x = numpy.array([data['spatial']['x'][i] for i in indices])
            numpy.random.shuffle(indices)
            y = numpy.array([data['spatial']['y'][i] for i in indices])

            data['spatial']['x'] = x
            data['spatial']['y'] = y
            print '#getData(): galaxy locations shuffled'


        kills = []
        killer = []

        min_colours = 2
        if overrides.has_key('min_colours'):
            min_colours = int(overrides['min_colours'])
        
        if len(limits.keys()) > 0 and (no_constraint == False and make_subs_bailsafe == False):
            if permit_dropouts:
#                print 'multi survey'
                for key in limits.keys():
                    no_detect = (data['magnitudes'][key] >= limits[key]) | (data['magnitudes'][key] < 0.0)
                    data['magnitudes'][key][no_detect] = numpy.nan
                
                filter = numpy.zeros(len(data['spatial']['x']))
                for col in original_colours:
                    m1_name = col.split('-')[0]
                    m2_name = col.split('-')[-1]
                    m1 = data['magnitudes'][col.split('-')[0]]
                    m2 = data['magnitudes'][col.split('-')[-1]]
                    col = ''.join(col.split('-'))
                    data['colours'][col] = m1-m2
                    f1 = numpy.isnan(data['colours'][col])
                    f1 = (f1 == False)
                    filter[f1] += 1
                    
                filter = (filter >= min_colours)
                constraint = constraint & filter
#                print 'filter done'

                    

#                filter = retain
#                constraint = constraint & filter

            else:
                for key in limits.keys():
                    filter = (data['magnitudes'][key] <= limits[key]) & (data['magnitudes'][key] > 0.0)
#                    test = (filter == False) & (constraint == True)
                    constraint = constraint & filter
#                    if len(kills) == 0:
#                        kills.append(countOf(test,True))
#                    else:
#                        kills.append(countOf(test,True))
#                        killer.append(key)
                        
        nmag_cut = countOf(constraint,True)
        if subset_index != None and sub_check and (grid != None or (selection != None and selection.has_key('grid'))):
            if grid == None:
                grid = selection['grid']


            if len(grid.keys()) == 1:
                if grid.keys()[0] in [-1,0]:
                    print '#getData(): Removing superflous outer shell in grid for a non-defined stripe'
                    grid = grid[grid.keys()[0]]

#            if stripe:
#                _grid = {subset_index:grid}
#                grid = _grid


            keys = grid.keys()
            sub_ids = [grid[id]['sub_id'] for id in keys]
            try:
                _this_grid = grid[keys[indexOf(sub_ids,subset_index)]]
            except:
                print '#getData(): Subset index not found in grid data - faking it!!'
                _this_grid = {}
                _this_grid['xmin'] = min(data['spatial']['x'])
                _this_grid['xmax'] = max(data['spatial']['x'])
                _this_grid['ymin'] = min(data['spatial']['y'])
                _this_grid['ymax'] = max(data['spatial']['y'])
                _this_grid['units'] = 'radians'
                _this_grid['sub_id'] = subset_index
                _this_grid['neighbourID'] = []
                xmin,xmax = _this_grid['xmin'],_this_grid['xmax']
                ymin,ymax = _this_grid['ymin'],_this_grid['ymax']
                xc = xmax-xmin
                yc = ymax-ymin
                n_dna = hash('%1.8f %1.8f'%(xc,yc))
                _this_grid['id'] = n_dna
                grid[subset_index] = _this_grid

#                grid = _this_grid
                
            xmin,xmax = _this_grid['xmin'],_this_grid['xmax']
            ymin,ymax = _this_grid['ymin'],_this_grid['ymax']

            if _this_grid['units'] == 'degrees':
                xmin = numpy.radians(xmin)
                xmax = numpy.radians(xmax)
                ymin = numpy.radians(ymin)
                ymax = numpy.radians(ymax)

            filter = (data['spatial']['x'] >= xmin) & (data['spatial']['x'] <= xmax) & (data['spatial']['y'] >= ymin) & (data['spatial']['y'] <= ymax)
            print '#getData(): Total # of input galaxies: %d, of which %d in survey RA limits'%(len(constraint),countOf(filter,True))

            n_spatial = countOf(filter,True)

            constraint = constraint & filter

 #           if stripe:
 #               grid = _grid.values()[0]
            

#            killer.append('spatial')
#            kills.append(countOf(filter,False))

        elif subset_index != None and sub_check and grid == None:
            from tools import which_limits
            xmin,xmax = which_limits(subset_index,wrap=wrap)
            filter = (data['spatial']['x'] >= xmin) & (data['spatial']['x'] <= xmax)
            print '#getData(): Defining subset limits from which_limits() - this is likely wrong..'
            print '#getData(): Total # of input galaxies: %d, of which %d in survey RA limits'%(len(constraint),countOf(filter,True))

            n_spatial = countOf(filter,True)

            constraint = constraint & filter
#            killer.append('spatial')
#            kills.append(countOf(filter,False))


        else:
            n_spatial = len(constraint)


        print '#getData(): %d input galaxies. %d are in survey RA limits. %d are in mag limits (%d both)'%(len(constraint),n_spatial,nmag_cut,countOf(constraint,True))

        if len(crange.keys()) == 0:
#            print 'Will do a fancy interactive thing here...'
            user_lims = {}
            crange = {}
#            for key in data['colours'].keys():
#                crange[key] = [-999,999]
#            string = instream[line_numbers['range']]
            
        elif len(crange.keys()) > 0 and not permit_dropouts and (no_constraint == False and make_subs_bailsafe == False):
            for key in crange.keys():
                _filter = (data['colours'][key] >= crange[key][0]) & (data['colours'][key] <= crange[key][1])
                constraint = constraint & _filter
                if (True in _filter) == False:
                    print '#getData(): Colour range %s knocked out selection'%(key)


        elif len(crange.keys()) > 0 and permit_dropouts and (no_constraint == False and make_subs_bailsafe == False):
            filter = numpy.zeros(len(data['spatial']['x']))
            for key in crange.keys():
#                f1 = (not numpy.isnan(data['colours'][key].all())) and (data['colours'][key] >= crange[key][0]) & (data['colours'][key] <= crange[key][1])
                f1 = (numpy.isnan(data['colours'][key]) == False) & (data['colours'][key] >= crange[key][0]) & (data['colours'][key] <= crange[key][1])
#                (data['colours'][key] >= crange[key][0]) & (data['colours'][key] <= crange[key][1])
                filter[f1] += 1
            filter = (filter >= min_colours)
            constraint = constraint & filter
#            print 'colour limit filter done'
                
            

#        if 1:
        if maxBCG_mode:
            from clusters import vertex,galaxy,vertex_list,point_inside_cell
            import numpy
            from tools import which_limits
            x = data['magnitudes']['i']
            y = data['colours']['gr']
            xmin,xmax = which_limits(-1)
            Ax,Ay = 13,1.67
            Bx,By = 13,0.79
            Cx,Cy = 17.3,0.79
            C2x,C2y = 18.85,1.1
            Dx,Dy = 19.8,1.42
            Ex,Ey = 19.8,1.67
            
            gal = galaxy(0,16,1.2,dummy=True)
            
            sel = {}
            sel['limits'] = {}
            sel['limits']['xmin'],sel['limits']['xmax'] = -10,21
            sel['limits']['ymin'],sel['limits']['ymax'] = 0.7,1.8
            
            vertices = [vertex(0,Ax,Ay),vertex(0,Bx,By),vertex(0,Cx,Cy),vertex(0,C2x,C2y),vertex(0,Dx,Dy),vertex(0,Ex,Ey)]
            selection_region = vertex_list(vertices,sel,gal,keep_all=True)

            constraint = numpy.array(numpy.ones(len(data['spatial']['x'])),dtype='bool')

            filter = (data['spatial']['x'] >= xmin) & (data['spatial']['x'] <= xmax)
            constraint = constraint & filter

            from dataIO import Progress
            vert_stat = Progress(len(x))
            keep = []
            print '#getData(): Determining which sources fall inside maxBCG selection'
            for i in xrange(len(x)):
                keep.append(point_inside_cell(x[i],y[i],selection_region.verts))
                vert_stat.update(0)
            
            keep = numpy.array(keep)

            constraint = constraint & keep

            update = True
            new_name = input_name.split('.fits')[0]+'_maxbcg'+'.fits'

            #point_inside_cell(x,y,selection_region.verts)
            

        if (update == True) and '_shuffle' not in input_name and (colour_shuffle or loc_shuffle) and '_shuffle' in new_name:
            constraint = numpy.array(numpy.ones(len(data['spatial']['x'])),dtype='bool')

        for key in data.keys():
            for sub_key in data[key]:
#                print key,sub_key
                if type(data[key][sub_key]) != type(None) and len(data[key][sub_key]) != 0:
#                if data[key][sub_key] != None and len(data[key][sub_key]) != 0:
                    data[key][sub_key] = data[key][sub_key][constraint]

                   

        if angular_mask:
            if not 'weight' in tbl.names:
                from tools import angular_mask
                angular_mask(data,plot=plot,frac=0.99)
                update = True
                new_name = input_name.split('.fits')[0]+'_extra'+'.fits'
                if plot:
                    try:
                        from dataIO import load
                        all = load('dxs_wt_024.dat')
                        from detections import plot_clusters
                        xlim = plt.xlim()
                        ylim = plt.ylim()
                        plot_clusters(all.clusters,rings=True,centroids=True,colour='!white',passback=True)
                        plt.xlim(xlim)
                        plt.ylim(ylim)
                    except:
                        pass
                    plt.show()
            else:
                data['extras']['weight'] = numpy.array(tbl.get_column('weight'))
            data['extras']['weight'] = data['extras']['weight'][constraint]
                
        for key in extinctions.keys():
            extinctions[key] = extinctions[key][constraint]
        if 1:
            for key in src_magnitudes.keys():
                src_magnitudes[key] = src_magnitudes[key][constraint]
        
        if ((subset_index != None) and (not os.path.isfile(new_name)) and sub_write) or update == True:
#            print data['spatial'].keys()
#            print data['extras'].keys()
#            raise SystemExit(0)

            import pyfits
            columns = []
            columns.append(pyfits.Column(name=ra_name,array=numpy.degrees(data['spatial']['x']),format='D'))
            columns.append(pyfits.Column(name=dec_name,array=numpy.degrees(data['spatial']['y']),format='D'))

            try:
                columns.append(pyfits.Column(name=sz_name,array=numpy.array(data['spatial']['specz']),format='D'))
            except:
                pass
            try:
                columns.append(pyfits.Column(name=pz_name,array=numpy.array(data['spatial']['photoz']),format='D'))
            except:
                pass

            for i in xrange(len(cid)): 
#                columns.append(pyfits.Column(name=cid[i],array=src_magnitudes[band_names[i]][constraint],format='D'))
                columns.append(pyfits.Column(name=cid[i],array=src_magnitudes[band_names[i]],format='D'))

            if len(eid) > 0:
                for i in xrange(len(eid)): 
#                    columns.append(pyfits.Column(name=eid[i],array=extinctions[eid[i]][constraint],format='D'))
                    columns.append(pyfits.Column(name=eid[i],array=extinctions[eid[i]],format='D'))

            if uid_name != '':
                columns.append(pyfits.Column(name=uid_name,array=data['extras']['id'],format='K'))

            _npflts = [numpy.float,numpy.float32,numpy.float64,numpy.float128]
            _npints = [numpy.int,numpy.int0,numpy.int8,numpy.int16,numpy.int32,numpy.int64]

            for col_entry in data['extras'].keys():
                if type(data['extras'][col_entry]) != type(None) and col_entry != 'id':
                    if countOf(data['extras'][col_entry],None) != len(data['extras'][col_entry]):
                        fmt = 'D'
                        if type(data['extras'][col_entry][0]) in _npflts or type(data['extras'][col_entry][0]) == type(1.1):
                            fmt = 'D'
                        elif type(data['extras'][col_entry][0]) in _npints or type(data['extras'][col_entry][0]) == type(1):
                            fmt = 'K'
                        elif type(data['extras'][col_entry][0]) == type('string'):
                            lmax = max([len(_s) for _s in data['extras'][col_entry]])
                            fmt = '%dA'%(lmax)
                            
                        columns.append(pyfits.Column(name=col_entry,array=data['extras'][col_entry],format=fmt))

            fits_defs = pyfits.ColDefs(columns)
            table = pyfits.new_table(fits_defs)
            print 'Writing file to %s'%(new_name)
            table.writeto(new_name,clobber=True)
            if colour_shuffle or loc_shuffle:
                print 3*'\n'
                print '#getData(): Shuffled. Change .cfg <filename> entry to point to %s'%(new_name)
                print 3*'\n'
                raise SystemExit(0)


    elif 'fits' in catalogue and dataset == 'mock':
#        import fits_table as fits
        tbl = fits.fits_table(catalogue)
        x,y = tbl.get_column('alpha'),tbl.get_column('delta')
        vx = tbl.get_column('vx')
        vy = tbl.get_column('vy')
        vz = tbl.get_column('vz')
        gmag = tbl.get_column('magpgo_tot_ex')
        rmag = tbl.get_column('magpro_tot_ex')
        imag = tbl.get_column('magpio_tot_ex')
        zmag = tbl.get_column('magpzo_tot_ex')
        ymag = tbl.get_column('magpyo_tot_ex')
        stellar_mass = tbl.get_column('mstars_tot')
        halo_mass = tbl.get_column('mhhalo')
        halo_id = tbl.get_column('nbody_halo_id')
        halo_idrep = tbl.get_column('idrep')
        halo_jm = tbl.get_column('jm')
        halo_jmtree = tbl.get_column('jtree')
        halo_zcos = tbl.get_column('z_obs')
        guess_limits = False


    data['x'] = data['spatial']['x']
    data['y'] = data['spatial']['y']
    data['shortcuts'] = {}
    data['shortcuts']['x'] = {}
    data['shortcuts']['y'] = {}




    #assume limits, colour list etc returned with data....


    start = (0.,0.,100.)
    finish = start
    change = (0.,0.,0.)
    combined = (start,finish,change)
    colours = (combined,combined,combined,combined)





    if args:
        dc = args['dc']
        probthresh = args['probthresh']
        rho_crit = args['rho_crit']
        smallest_cluster = args['smallest_cluster']
        halomass = args['halomass']
        metric = args['metric']
        outlier_metric = args['outlier_metric']
        outlier_no_clobber = args['outlier_no_clobber']
        
        if len(crange.keys()) > 0:
            range = crange
        else:
            range = {'gr':[0.47,2.00+dc/2.0,dc],'ri':[0.0,1.5+dc/2.0,dc],'iz':[-0.1,1.1+dc/2.0,dc]}

        x,y = data['spatial']['x'],data['spatial']['y']


        empty_cell = False

        if dataset == 'swire' or dataset == 'mds' or dataset == 'hershel':
            args['win'] = max((max(x)-min(x)),(max(y)-min(y)))/2.0
        else:
            try:
                args['win'] = min((max(x)-min(x)),(max(y)-min(y)))/2.0
            except ValueError:
                pass                                                 ##x or y are empty
                
#        phot_err = 27.5

        core = classes.detections(None,None,dummy=True)
        core.selection = {}


#        selection = core.selection
        core.selection['limits'] = {}
        try:
            core.selection['limits']['xmin'],core.selection['limits']['xmax'] = min(x),max(x)
            core.selection['limits']['ymin'],core.selection['limits']['ymax'] = min(y),max(y)
        except ValueError:
            print 'I think this might be an empty cell...'
            empty_cell = True
            core.selection['limits']['xmin'],core.selection['limits']['xmax'] = 1,2
            core.selection['limits']['ymin'],core.selection['limits']['ymax'] = 1,2
            pass
            
        core.selection['limits']['colours'] = {}
        core.selection['limits']['magnitudes'] = {}
        core.selection['limits']['cmr_range'] = range
        for key in range:
            core.selection['limits']['colours'][key] = range[key]


        for band in band_names:
            if limits[band] != 40.0:
                core.selection['limits']['magnitudes'][band] = {'limit':limits[band]}
            else:
                core.selection['limits']['magnitudes'][band] = {'limit':99}

        core.selection['io'] = {}
        core.selection['io']['dataset'] = dataset
        core.selection['io']['sourcefile'] = catalogue
        core.selection['io']['permit_dropouts'] = permit_dropouts
        if input_name != catalogue:
            core.selection['io']['sourcefile'] = input_name
        #now do md5 hash on the source data?
        if md5:
            from tools import md5file
            print 'Calculating MD5 hash for sourcefile'
            core.selection['io']['md5hash'] = md5file(catalogue)
            #I think we also want the args soecified at terminal as well
            import sys
            core.selection['io']['cmd'] = ' '.join(sys.argv)
            
        core.selection['io']['cfgsrc'] = cfg
        core.selection['io']['observed_data'] = not 'mock' in dataset
        core.selection['io']['wrap'] = wrap
        
        

        try:
            from bzrlib.branch import BzrBranch
            branch =  BzrBranch.open('.')
            orca_version = branch.last_revision_info()[0]
            core.selection['io']['orca_version'] = orca_version
        except:
            pass
        if data['extras'].has_key('offset') and max(data['extras']['offset']) > 0:
            core.selection['io']['wrap'] = True
            if offset_threshold:
                core.selection['io']['wrap_threshold'] = offset_threshold
        
        if grid != None:
            local_grid = {}
            keys = grid.keys()
#            print grid.keys()
#            print grid[keys[0]]
            nested = False
            if len(keys) == 1 and not 'sub_id' in grid[keys[0]].keys():
                sub_ids = [grid[keys[0]][id]['sub_id'] for id in grid[keys[0]].keys()]
                nested = True
                keys = grid[keys[0]].keys()
            else:
                sub_ids = [grid[id]['sub_id'] for id in keys]
                
            try:
                try:
                    if nested:
                        _grid = grid.values()[0][keys[indexOf(sub_ids,subset_index)]]
                    else:
                        _grid = grid[keys[indexOf(sub_ids,subset_index)]]
                except ValueError:
                    if 'self' in keys:
                        _grid = grid['self']
                neighbour_id = _grid['neighbourID']
                local_grid['self'] = _grid
                for n in neighbour_id:
                    local_grid[grid[n]['sub_id']] = grid[n]

                core.selection['grid'] = local_grid
            except ValueError:
                pass


        core.selection['detection'] = {'probthresh':probthresh,'rho_crit':rho_crit,'smallest_cluster':smallest_cluster,'metric':metric,'subset_index':subset_index,'sequence':sequence}
        core.selection['detection']['detection_area'] = element_area(degrees(core.selection['limits']['xmin']),degrees(core.selection['limits']['xmax']),degrees(core.selection['limits']['ymin']),degrees(core.selection['limits']['ymax']))
#        core.selection['detection']['detection_area'] = (degrees(core.selection['limits']['xmax']-core.selection['limits']['xmin'])*degrees(core.selection['limits']['ymax']-core.selection['limits']['ymin']))}
        core.selection['detection']['mean_density'] = 1.0/core.selection['detection']['detection_area']
#        print selection['limits']['xmin'],selection['limits']['xmax']
#        print selection['limits']['ymin'],selection['limits']['ymax']
#        print selection['detection']['detection_area']
#        raise SystemExit(0)
        
        if subset_index == None:
            core.selection['detection']['subset_index'] = None
        if stripe:
            core.selection['detection']['stripe'] = stripe

        core.selection['overrides'] = {}
        if len(overrides.keys()) > 0:
            for okey in overrides.keys():
                core.selection['overrides'][okey] = overrides[okey]
        core.selection['colours'] = {}
        for col in data['colours']:
            core.selection['colours'][col] = {'intercept':0.0,'slope':0.0,'width':100.0}
            core.selection['colours'][col]['fitfuncs'] = {}
            for _name in band_names:
                if col.split(_name)[-1] == '':
                    core.selection['colours'][col]['magnitude'] = _name

        if len(filter_slopes.keys()) > 0:
            for col in filter_slopes.keys():
                core.selection['colours'][col]['slope'] = filter_slopes[col]

        if len(filter_widths.keys()) > 0:
            for col in filter_widths.keys():
                core.selection['colours'][col]['width'] = filter_widths[col]




        junk = """clean this up!!!! selection? core.selection? Shouldn't this be the same, and then returned back? FFS"""



        core.selection['mock_data'] = {'halomass':halomass,'outlier_metric':outlier_metric,'no_clobber':outlier_no_clobber}
        

        core.selection['defaults'] = {'colour':core.selection['detection']['sequence'][0]}
        col = core.selection['defaults']['colour']
        core.selection['defaults']['magnitude'] = core.selection['colours'][col]['magnitude']
        mag = core.selection['defaults']['magnitude']
        if phot_err and core.selection['limits']['magnitudes'].has_key(mag):
            core.selection['limits']['magnitudes'][mag]['photometric_error'] = phot_err
        

        for key in data['magnitudes'].keys():
            if 1:
                _mag = key
                try:
                    mag_lim = core.selection['limits']['magnitudes'][_mag]['limit']
                except:
                    mags = deepcopy(data['magnitudes'][_mag])
                    filter = (mags < 50) & (mags > 0)
                    mags = mags[filter]
                    mags = list(mags)

                    mags.sort()
                    mags.reverse()
                    faint_lim = median(mags[:min(250,len(mags))])
                    mags.reverse()
                    bright_lim = median(mags[:min(250,len(mags))])
                    del mags
                    core.selection['limits']['magnitudes']['!'+_mag] = {}
                    core.selection['limits']['magnitudes']['!'+_mag]['limit'] = faint_lim
        for key in sequence:
            if not key in core.selection['limits']['cmr_range']:
                _mag = core.selection['colours'][key]['magnitude']
                cols = deepcopy(data['colours'][key])
                mags = deepcopy(data['magnitudes'][_mag])

                try:
                    maglim = core.selection['limits']['magnitudes'][_mag]['limit']
                except:
                    maglim = core.selection['limits']['magnitudes']['!'+_mag]['limit']

                filter = (cols < 5.0) & (cols > -5.0) & (mags < maglim -1) & (mags > maglim -2)

                cols = cols[filter]
                m,s = median(cols),std(cols)
#                core.selection['limits']['cmr_range']['!'+key] = [m-(2*s),m+(3*s),dc]
                core.selection['limits']['cmr_range'][key] = [m-(2*s),m+(3*s),dc]
                
        if core.selection.has_key('grid'):
            core.set_limits()


                    
        if plot == False:
            return data,colours,tbl,range,core            
#            return data,colours,None,range,core

    if plot:
        guess_limits = False
        from tools import rows_cols

        fig = plt.figure()
        colours = core.selection['detection']['sequence']

        plot = iter(arange(len(colours))+1)
        n = len(colours)
        nrows,ncols = rows_cols(n)

        plots = iter(arange(len(colours))+1)

        for c in colours:
            mag = core.selection['colours'][c]['magnitude']
            try:
                xmax = core.selection['limits']['magnitudes'][mag]['limit']
            except:
                xmax = 24.0
            try:
                xmin = 13
                ymin = core.selection['limits']['cmr_range'][c][0]
                ymax = core.selection['limits']['cmr_range'][c][1]

            except:
                xmin = 13
                ymin = core.selection['limits']['cmr_range']['!'+c][0]
                ymax = core.selection['limits']['cmr_range']['!'+c][1]
                

#            except:
#                print 'Guessing limits for colour %s'%(c)
#                xmin = 13
#                ymin = 0.0
#                ymax = 2.0
                


            _col = data['colours'][c]
            _mag = data['magnitudes'][mag]
            constraint = (_col <= ymax) & (_col >= ymin) & (_mag <= xmax) & (_mag >= xmin)
            _col = _col[constraint]
            _mag = _mag[constraint]

            plt.subplot(nrows,ncols,plot.next())
            plt.hexbin(_mag,_col,gridsize=50,bins='log',cmap=plt.cm.jet,antialiased=True)
            plt.xlabel(mag)
            plt.ylabel(c)
            plt.xlim(xmin,xmax)
            plt.ylim(ymin,ymax)

        plt.subplots_adjust(wspace=0.3,right=0.97,left=0.10,bottom=0.07,top=0.97,hspace=0.3)
        plt.show()
        if args:
            return data,colours,tbl,range,core
            
    return data
