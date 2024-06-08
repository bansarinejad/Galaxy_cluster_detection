from operator import countOf,indexOf
from clusters import galaxy,cluster
from numpy import compress,equal,hypot,degrees,radians,log10,arange,unique,median,average,hypot,mean,compress,zeros,sort,poly1d,sqrt,std,array,less,exp
from copy import copy,deepcopy
import dataIO as datacuts
import gc
from science import science
#from tools import common_area,objectify,all_good_detected,count_names
import os
#if os.getenv('HOSTNAME') != 'dubris':
#if os.getenv('HOSTNAME') == None or os.getenv('HOSTNAME') != 'dubris' and not (os.getenv('HOSTNAME') in ['pollux','castor','m3001','m3002'] or os.getenv('HOSTNAME')[0] == 'm'):
if os.getenv('HOSTNAME') == None or os.getenv('HOSTNAME') != 'dubris' and not (os.getenv('HOSTNAME') in ['m3001','m3002'] or os.getenv('HOSTNAME')[0] == 'm'):
    import matplotlib.pyplot as plt
    import matplotlib
from parallel_detect import parallel_detect
from interaction import interaction
from math import pi
import numpy
from numpy import radians,degrees,std
from parallel_qhull import qhull
from parallel_vertices import deploy_vertices


class detections(object):
    limits = {'ri':{'xlim':[15.5,24],'ylim':[-0.2,2.0]},'gr':{'xlim':[15.5,24],'ylim':[-0.2,2.3]},\
                  'iz':{'xlim':[15.5,24],'ylim':[-0.2,1.2]},'zy':{'xlim':[15.5,24],'ylim':[-0.2,0.6]}}


#class detections(science):


    def __init__(self,cuts,selection,make_perfect=False,sourcedata=None,cherrypick=False,set_density=False,dummy=False):
        #Metrics
        self.merge_crits = {1:0,2:0,3:0,4:0,5:0}
        self.comp = 0.
        self.pure = 0.
        self.qual = 0.
        self.acc_r = 0.
        self.acc_sm = 0.
        self.acc_v = 0.
        self.frag = 0.
        self.wrap_mode = 'detected'
        self.prev_wrap_mode = None

        #Data from co-detection in two colours
        self.colour_colour = None
        self.co_detections = False

        self.source_data = None
        self.source_galaxies = cuts
        self.cuts = cuts
        self.mia = []
        # these are cluster names not found in the perfect list...
        self.cluster_galaxies = None
        self.perfect_cluster_galaxies = []
        self.clusters = None
        self.perfect_clusters = []
        self.field_galaxies = []


        #catalogue of external data to plot
        self.external_data = None

        # can this detection be used to pick the best clusters from?
        self.cherrypick = cherrypick
        # flavour: default is cherrypicking among other detection classes
        # within the same colour space. This variable is overridden by
        # cherrypick from 'intra' to 'inter' when cherrypicking from
        # clusters in different colour spaces (eg g-r and r-i cluster detections)
        self.cherry_flavour = ''
        if self.cherrypick: self.cherry_flavour = 'intra'
        self.cherry_flavour = 'intra'

        if dummy and not selection:
            return


        self.selection = deepcopy(selection)
        try:
            self.sub_id = self.selection['detection']['subset_index']
        except KeyError:
            self.sub_id = None


        if selection.has_key('overrides') and selection['overrides'].has_key('area'):
            from math import pi
            self.selection['detection']['detection_area_radians'] = selection['overrides']['area']
            self.selection['detection']['detection_area'] = selection['overrides']['area']*((360.0*360.0)/(4*pi*pi))

        else:
            self.selection['detection']['detection_area_radians'] = float((self.selection['limits']['xmax']-self.selection['limits']['xmin'])*(self.selection['limits']['ymax']-self.selection['limits']['ymin']))
            self.selection['detection']['detection_area'] = (degrees(self.selection['limits']['xmax']-self.selection['limits']['xmin'])*degrees(self.selection['limits']['ymax']\
                                                                                                           -self.selection['limits']['ymin']))

            try:
                self.selection['detection']['detection_area'] = selection['detection']['detection_area']
            except:
                pass

            try:
                self.selection['detection']['detection_area_radians'] = selection['detection']['detection_area_radians']
            except KeyError:
                _a = self.selection['detection']['detection_area']
                from math import pi
                self.selection['detection']['detection_area_radians'] = 4*pi*pi*_a/(360.0*360.0)
                

                
#        self.selection['detection']['probthresh_function'] = poly1d([ -4.94880000e-22,   2.53161646e-14,  -4.28928452e-07,  -2.23225549e+01])
#        self.selection['probthresh'] = self.selection['probthresh_function'](self.selection['mean_density_est'])
#        self.selection['detection']['rc_fits'] = [1.71034706e+05,6.74099161e-07,-3.43812816e+05,2.13349188e+00]

#        self.rc_fit = self.selection['detection']['rc_fits']
#        from math import exp
#        lambda p, x: p[3]+x/p[0]*exp(-p[1]*(x+p[2]))

        if dummy:
            return

        self.selection['detection']['ngal'] = len(self.cuts['x'])
        self.selection['detection']['mean_density'] = float(self.selection['detection']['ngal'])/float((self.selection['limits']['xmax']-self.selection['limits']['xmin'])*(self.selection['limits']['ymax']-self.selection['limits']['ymin']))

    
    def __getitem__(self,index):
        return self.clusters[index]


    def diff(self,det2):
        """
        Do a comparison of this detection with det2
        """
        ignores = []
        ignores.append('colours|*|intercept')
        ignores.append('colours|*|slope')
        ignores.append('colours|*|width')
        ignores.append('colours|*|trainer_width')

        
        base = self
        other = det2

        base_sel = base.selection
        other_sel = other.selection

        
        
        #0th level: same keys?
        for key in base.selection.keys():
            if not key in other.selection.keys():
                print 'Missing in OTHER.selection: %s'%(key)
                ignores.append(key)
        for key in other.selection.keys():
            if not key in base.selection.keys():
                print 'Missing in BASE.selection: %s'%(key)
                ignores.append(key)

        #1st level: check dictionary values against each other

        checked = []

#        for key in base.selection.keys()


        
        
    def publish(self,url='auto',clobber=False,grid=False,allow_blank=False,bcg_centred=True,pan_plot=True,hdf5=True,fits=True,uber=False):
        import os,shutil
        self.wrap_toggle(mode='observed')
#        url = 'http://astro.dur.ac.uk/~dmurphy/orca/data/panstarrs/md08/images/'
        
        if uber:
            source_clusters = self.clusters
            try:
                ucg = self.uber_cluster_galaxies[0]
            except:
                self.retrieve_associates(uber=True,rfrac=0.8)
            self.clusters = self.uber_clusters

        img_srcs = [cl.img_name(bcg_mode=bcg_centred) for cl in self.clusters]
        if grid:
            img_gsrc = []
            for img in img_srcs:
                img_gsrc.append('/'.join(img.split('/')[:-1])+'/grid/'+img.split('/')[-1])
            img_srcs = img_gsrc

                
        img_test = [os.path.isfile(image) for image in img_srcs]
        incomplete = False
        if False in img_test and (not allow_blank):
            gg = ''
            if grid:
                gg = 'grid '

            print 'At least one of the clusters in this class has no %simage - run self.retrieve_images(no_clobber=True) first'%(gg)
        
        elif False in img_test:
            print 'Warning - not all clusters have images'
            incomplete = True



        if url == 'auto':
            url = 'http://astro.dur.ac.uk/~dmurphy/orca/%s/'%(self.selection['io']['dataset'])
        elif url[-1] != '/':
            url = url+'/'


        storage = './%s/'%(self.selection['io']['dataset'])
        if os.path.isdir(storage) and not clobber:
            print 'Removing directory'
            os.system('rm -Rf %s'%(storage))
#            os.rmdir(img_dir)
        os.mkdir(storage)
        
        img_dir = storage+'images/'
        if os.path.isdir(img_dir) and not clobber:
            print 'Removing directory'
            os.system('rm -Rf %s'%(img_dir))
#            os.rmdir(img_dir)
        os.mkdir(img_dir)

        if pan_plot:
            pan_store= storage+'membership/'
            if os.path.isdir(pan_store) and not clobber:
                print 'Removing directory'
                os.system('rm -Rf %s'%(pan_store))
            os.mkdir(pan_store)

        for i in xrange(len(self.clusters)):
            cl = self.clusters[i]
            if incomplete:
                if not os.path.isfile(cl.img_name(bcg_mode=bcg_centred)):
                    cl.image = ''
                    continue
            if pan_plot:
                if incomplete and not os.path.isfile(cl.img_name(bcg_mode=bcg_centred)):
                    pass
                else:
                    cl.pan_plot(loc=pan_store,store=True,bcg_centred=bcg_centred,assoc=True)
                    print 'Panplot generated'
            src = img_srcs[i]
            dest = img_dir+'%d.jpg'%(abs(cl.dna))
            shutil.copy(src,dest)
            print 'Image to %s'%(dest)
            cl.image = url+'images/%d.jpg'%(abs(cl.dna))
        print 'Cluster images copied, cl.image strings updated'


        if hdf5:
            outname = storage+'%s.hdf5'%(self.selection['io']['dataset'])
            self.to_hdf5(output=outname,toplevel_fits=fits)
        if uber:
            self.clusters = source_clusters

        self.wrap_toggle()

    def to_kml(self,output='out.kml',mode='clusters'):
        self.wrap_toggle(mode='observed')
        import numpy
        header = """<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom" hint="target=sky">
<Document>
	<name>out.kml</name>
	<open>1</open>
	<Style id="clusterBCG">
		<BalloonStyle>

<text><![CDATA[
<center>Brightest member of cluster: <b><br>$[name]</br></b></center><br></br>
<table width="60%" style="background-color:white" border="1" cellpadding="2" cellspacing="0" align = "center"> 
	<tr> 
		<td>Galaxy#</td> 
		<td>$[posn]</td> 
	</tr> 
	<tr> 
		<td>specz</td> 
		<td>$[zspec]</td> 
	</tr> 
	<tr> 
		<td>photoz</td> 
		<td>$[zphoto]</td> 
	</tr> 
</table> 

]]>
</text>

		</BalloonStyle>
	</Style>
	<Style id="sn_red-pushpin">
		<IconStyle>
			<color>66ffffff</color>
			<scale>1.1</scale>
			<Icon>
				<href>http://maps.google.com/mapfiles/kml/pushpin/red-pushpin.png</href>
			</Icon>
			<hotSpot x="20" y="2" xunits="pixels" yunits="pixels"/>
		</IconStyle>
		<ListStyle>
		</ListStyle>
<BalloonStyle>
<text><![CDATA[
<center>Brightest member of cluster: <b><br>$[name]</br></b></center><br></br>
<table width="60%" style="background-color:white" border="1" cellpadding="2" cellspacing="0" align = "center"> 
	<tr> 
		<td>Galaxy#</td> 
		<td>$[posn]</td> 
	</tr> 
	<tr> 
		<td>specz</td> 
		<td>$[zspec]</td> 
	</tr> 
	<tr> 
		<td>photoz</td> 
		<td>$[zphoto]</td> 
	</tr> 
</table> 

]]>
</text>

		</BalloonStyle>


	</Style>
	<Style id="clusterBCG1">
		<IconStyle>
			<color>6600ffff</color>
			<scale>2</scale>
		</IconStyle>
		<LabelStyle>
			<scale>1.5</scale>
		</LabelStyle>
		<BalloonStyle>
<text><![CDATA[
<center>Member of cluster: <b><br>$[name]</br></b></center><br></br>
<table width="60%" style="background-color:white" border="1" cellpadding="2" cellspacing="0" align = "center"> 
	<tr> 
		<td>Galaxy#</td> 
		<td>$[posn]</td> 
	</tr> 
	<tr> 
		<td>specz</td> 
		<td>$[zspec]</td> 
	</tr> 
	<tr> 
		<td>photoz</td> 
		<td>$[zphoto]</td> 
	</tr> 
</table> 

]]>
</text>
		</BalloonStyle>
	</Style>
	<Style id="clusterBCG2">
		<IconStyle>
			<color>6600ffff</color>
			<scale>2</scale>
		</IconStyle>
		<LabelStyle>
			<scale>1.5</scale>
		</LabelStyle>
		<BalloonStyle>

<text><![CDATA[
<center>Member of cluster: <b><br>$[name]</br></b></center><br></br>
<table width="60%" style="background-color:white" border="1" cellpadding="2" cellspacing="0" align = "center"> 
	<tr> 
		<td>Galaxy#</td> 
		<td>$[posn]</td> 
	</tr> 
	<tr> 
		<td>specz</td> 
		<td>$[zspec]</td> 
	</tr> 
	<tr> 
		<td>photoz</td> 
		<td>$[zphoto]</td> 
	</tr> 
</table> 

]]>
</text>
		</BalloonStyle>
	</Style>
	<StyleMap id="clusterMember">
		<Pair>
			<key>normal</key>
			<styleUrl>#sn_red-pushpin</styleUrl>
		</Pair>
		<Pair>
			<key>highlight</key>
			<styleUrl>#sh_red-pushpin</styleUrl>
		</Pair>
	</StyleMap>
	<Style id="sh_red-pushpin">
		<IconStyle>
			<color>66ffffff</color>
			<scale>1.3</scale>
			<Icon>
				<href>http://maps.google.com/mapfiles/kml/pushpin/red-pushpin.png</href>
			</Icon>
			<hotSpot x="20" y="2" xunits="pixels" yunits="pixels"/>
		</IconStyle>
		<ListStyle>
		</ListStyle>
<BalloonStyle>
<text><![CDATA[
<center>Member of cluster: <b><br>$[name]</br></b></center><br></br>
<table width="60%" style="background-color:white" border="1" cellpadding="2" cellspacing="0" align = "center"> 
	<tr> 
		<td>Galaxy#</td> 
		<td>$[posn]</td> 
	</tr> 
	<tr> 
		<td>specz</td> 
		<td>$[zspec]</td> 
	</tr> 
	<tr> 
		<td>photoz</td> 
		<td>$[zphoto]</td> 
	</tr> 
</table> 

]]>
</text>
		</BalloonStyle>
	</Style>
	<StyleMap id="clusterBCG0">
		<Pair>
			<key>normal</key>
			<styleUrl>#clusterBCG1</styleUrl>
		</Pair>
		<Pair>
			<key>highlight</key>
			<styleUrl>#clusterBCG2</styleUrl>
		</Pair>
	</StyleMap>

"""


        footer = """</Document>
</kml> 
"""

#        print header
        def placemark(cluster,mode='cluster'):
            import numpy
            from math import sin,cos,degrees
            cl = cluster
            x,y = cluster.x,cluster.y
            r80 = 2.0*degrees(cl.extent)*60.0*60.0
            longitude = degrees(cl.bcg.x)-180.0
            latitude = degrees(cl.bcg.y)
            R = 6.378e6
            k = 1.1917536
            LookAtRange = (R*(k*sin(r80/2.0) - cos(r80/2.0) + 1))/1000.0
            gal = cl.bcg
            posn = str(indexOf(cluster.galaxies,gal))
            if not numpy.isnan(gal.extras['specz']):
                zspec = str('%1.3f'%(gal.extras['specz']))
            else:
                zspec = '--'
                
            if not numpy.isnan(gal.extras['photoz']):
                photoz = str('%1.3f'%(gal.extras['photoz']))
            else:
                photoz = '--'
            
            kml = """<Placemark>
<name>%s</name>
<description>
<![CDATA[
     This is the cluster BCG from PanSTARRS
      ]]>
    </description>
    <LookAt>
      <longitude>%f</longitude>
      <latitude>%f</latitude>
      <altitude>0</altitude>
      <range>%f</range>
      <tilt>0</tilt>
      <heading>0</heading>
    </LookAt>
    <styleUrl>#clusterBCG</styleUrl>
<ExtendedData>
        <Data name="zspec">
          <displayName><![CDATA[
            <b><i>The yardage is </i></b>
          ]]></displayName>
          <value>%s</value>
        </Data>
       <Data name="zphoto">
          <displayName><![CDATA[
            <b><i>The yardage is </i></b>
          ]]></displayName>
          <value>%s</value>
        </Data>
       <Data name="posn">
          <displayName><![CDATA[
            <b><i>The yardage is </i></b>
          ]]></displayName>
          <value>%s</value>
        </Data>
      </ExtendedData>
    <Point>
      <coordinates>%f,%f,0</coordinates>
   </Point>
</Placemark>
            """%(cluster.iau_name,longitude,latitude,LookAtRange,zspec,photoz,posn,longitude,latitude)



            points = ''
            import numpy
            print mode
            for gal in cluster.galaxies:
                if gal != cluster.bcg and mode != 'clusters':
                    longitude = degrees(gal.x)-180.0
                    latitude = degrees(gal.y)
                    posn = str(indexOf(cluster.galaxies,gal))
                    if not numpy.isnan(gal.extras['specz']):
                        zspec = str('%1.3f'%(gal.extras['specz']))
                    else:
                        zspec = '--'

                    if not numpy.isnan(gal.extras['photoz']):
                        photoz = str('%1.3f'%(gal.extras['photoz']))
                    else:
                        photoz = '--'

                    points = points + """<Placemark>
<name>%s</name>
<description>
<![CDATA[
     This is a cluster member from PanSTARRS
      ]]>
    </description>
    <LookAt>
      <longitude>%f</longitude>
      <latitude>%f</latitude>
      <altitude>0</altitude>
      <range>%f</range>
      <tilt>0</tilt>
      <heading>0</heading>
    </LookAt>
    <styleUrl>#clusterMember</styleUrl>
<ExtendedData>
        <Data name="zspec">
          <displayName><![CDATA[
            <b><i>The yardage is </i></b>
          ]]></displayName>
          <value>%s</value>
        </Data>
       <Data name="zphoto">
          <displayName><![CDATA[
            <b><i>The yardage is </i></b>
          ]]></displayName>
          <value>%s</value>
        </Data>
       <Data name="posn">
          <displayName><![CDATA[
            <b><i>The yardage is </i></b>
          ]]></displayName>
          <value>%s</value>
        </Data>
      </ExtendedData>
    <Point>
      <coordinates>%f,%f,0</coordinates>
   </Point>
</Placemark>
 """%(cluster.iau_name,longitude,latitude,LookAtRange,zspec,photoz,posn,longitude,latitude)
            kml = kml + points
            return kml


        for cl in self.clusters:
            header = header + placemark(cl,mode=mode)
        kml = header + footer
        file = open(output,'w')
        file.writelines(kml)
        file.close()
        self.wrap_toggle()
        return kml

    

    def to_vstcutouts(self,filter='gri',output='orca_vstcutouts.dat',centre='bcg',location='./orca_cutouts'):
        """\
        det.to_vstcutouts(self, => void

        \t Write ORCA cluster positions to a VSTcutout file to be processed

        """
        if location[-1] != '/':
            location += '/'
        stream_out = []
        if type(filter) == type([]):
            filter = ','.join(filter)
        
            
        from math import degrees
        for cluster in self.clusters:
            src = cluster
            if centre == 'bcg':
                src = cluster.bcg
            ra = degrees(src.x)
            dec = degrees(src.y)
            size = int(2*degrees(cluster.max_radius(with_galaxies=True))*60.0*60.0)
            name = location+str(cluster.dna)
            stream_out.append('%f %f %d %s %s\n'%(ra,dec,size,name,filter))
        file = open(output,'w')
        file.writelines(stream_out)
        file.close()
        print 'VST cutout file written to %s'%(output)
            

            
            

    def to_aladin(self,id=None):
        try:
            from astropy.vo.samp import SAMPIntegratedClient
        except:
            print 'Could not import from Astropy - is it installed?'
            return
        client = SAMPIntegratedClient()
        client.connect()
        client_list = client.get_registered_clients()
        aladin_sessions = []
        aladin_ids = []
        for cl_id in client_list:
            cl_meta = client.get_metadata(cl_id)
            if cl_meta.has_key('samp.name') and cl_meta['samp.name'] == 'Aladin':
                info = {'aladin.version':cl_meta['aladin.version'],'id':cl_id}
                aladin_sessions.append(info)
                aladin_ids.append(cl_id)
        if ((len(aladin_ids) > 1) and (id == None)) or ((len(aladin_ids) > 1) and (id != None) and (id not in aladin_ids)):
            if id != None:
                print """Specified id (%s) is not defined as an Aladin client - please specify one in the args (det.to_aladin(id='c3')):"""
            else:
                print """More than one instance of Aladin - in future, please specify one in the args (det.to_aladin(id='c3')):"""
            for cl in aladin_sessions:
                print 'ID = %s \t Aladin %s'%(cl['id'],cl['aladin.version'])

            from dataIO import question
            string = "Please enter the ID of the Aladin instance you want the ORCA detection transmitted to: "
            answer_type = 'string'
            mixed_mode = True
            escape = "Q"
            user_id = ''
            while user_id not in aladin_ids:
                user_id = question(string,answer_type,escape=escape,mixed_mode=mixed_mode)
            id = user_id

            
        else:
            id = aladin_ids[0]
                
        
        #write file to /tmp/something.fits
        #push this file to Aladin
        self.to_fits(output='/tmp/out.fits')
        fileURL = 'file:///tmp/out_galaxies.fits'
        params = {}
        params["url"] = fileURL
        params["name"] = "ORCA catalogue"
        message = {}
        message["samp.mtype"] = "table.load.fits"
        message["samp.params"] = params
        client.notify(id, message)
        print "Done"

    def to_fits(self,output='out.fits',img_src=False,custom=None,thesis=False,prefix='ORC',texify=False,extras=[],galcat=True,gals_only=False):
        """\
        det.to_fits(self,output='out.fits',img_src=False,custom=None,thesis=False,prefix='ORC',texify=False,extras=[],galcat=True,gals_only=False) => void

        \t Write ORCA cluster detection object to fits file(s).
        
        output     :  filename required for output
        img_src    :  include the [local] URL to the imaging data for this cluster (if exists)*
        custom     :  overrides the cluster entries entirely to those specified
        thesis     :  sets table parameters for a document written a long time ago in a galaxy far far away...
        prefix     :  IAU prefix requested
        texify     :  create a .tex table output
        extras     :  add extra keywords here (for cluster component only)
        galcat     :  create a fits file with galaxy data as well (cross-match with clusterID)
        gals_only  :  only output the galaxy catalogue


        * this field can be used in combination with the Topcat "Activation Action" to view stored cluster images
        """



        from tools import cat
        from copy import deepcopy
        self.wrap_toggle(mode='observed')
        cl = self.clusters[0]
        cl.parent = None
        cl0 = deepcopy(cl)
        cl.parent = self


        if galcat:
            detection = self    ##lazy - some code copied from tools.py:cuts()
            if (gals_only):
                outname = output
            else:
                outname = output
                if not '.fit' in outname:
                    print '#detection.to_fits(): warning - file output does not have .fits extension (will continue)'
                    ending = outname.split('.')[-1]
                    outname = '.'.join(outname.split('.')[:-1])+'_galaxies.%s'%(ending)
                    print 'This is the galaxy file: %s'%(outname)
                else:
                    outname = output.split('.fit')[0]+'_galaxies.fits'
            galaxies = self.cluster_galaxies
            mags = self.selection['limits']['magnitudes'].keys()
            galcat_extras=['specz','photoz']
            if len(detection.associated_cluster_galaxies) > 0:
                galcat_extras.append('assoc')
                for g in galaxies:
                    g.assoc = False
                for g in list(detection.associated_cluster_galaxies):
                    g.assoc = True
                galaxies = list(galaxies)+list(detection.associated_cluster_galaxies)

            if 1:
#            try:
                is_bcg = detection.cluster_galaxies[0].is_bcg
                galcat_extras.append('is_bcg')
            else:
#            except:
                a = 1

            from tools import gallist_cat
#            print outname
#            print outname
            out = gallist_cat(galaxies,output=outname,extras=galcat_extras,mags=mags,name_override=True,label=None)
            if gals_only:
                print 'Output only galaxies'
                return




        preserve = ['scatter']

        if not custom:
            custom = ['ID','ra','dec','ra_bcg','dec_bcg','ngal','redshift','redshift_code','theta80','concentration']
            custom = ['iau_name','RA','DEC','redshift','redshift_code','N_gal','bgc','scatter','theta80','concentration']



        if thesis:
            custom = ['iau_name','RA','DEC','cluster_z','cz_type','N_gal','bgc','scatter','theta80','concentration']
            for cl in self.clusters:
                cl.get_iau_name(prefix=prefix)
                if cl.x < 0:
                    cl.RA = degrees(cl.x) + 360
                else:
                    cl.RA = degrees(cl.x)
                cl.DEC = degrees(cl.y)
                cl.cluster_z = cl.z
                cl.cz_type = cl.z_type
                cl.N_gal = cl.ngalaxies
                try:
                    cl.bgc = cl.rich_data['bgc']['bgc']
                    cl.bgc = cl.agc1i
                except:
                    cl.bgc = numpy.nan
                cl.scatter = cl.cmr(fit=True,colour=cl.basis)[1]
                cl.theta80 = cl.r_x(frac=0.8)
                cl.concentration = cl.r_x(frac=0.8)/cl.r_x(frac=0.2)
        else:
            for cl in self.clusters:
                cl.r_method = cl.redshift

                if cl.bcg.x < 0:
                    cl.ra_bcg = degrees(cl.bcg.x) + 360
                else:
                    cl.ra_bcg = degrees(cl.bcg.x)
                cl.dec_bcg = degrees(cl.bcg.y)

                if 'bgc' in custom:
                    try:
                        cl.bgc = cl.rich_data['bgc']['bgc']
                    except:
                        cl.bgc = numpy.nan
                if 'ID' in custom:
                    cl.ID = abs(cl.dna)
                if 'RA' in custom:
                    if cl.x < 0:
                        cl.RA = degrees(cl.x) + 360
                    else:
                        cl.RA = degrees(cl.x)
                elif 'ra' in custom:
                    if cl.x < 0:
                        cl.ra = degrees(cl.x) + 360
                    else:
                        cl.ra = degrees(cl.x)
                if 'DEC' in custom:
                    cl.DEC = degrees(cl.y)
                elif 'dec' in custom:
                    cl.dec = degrees(cl.y)
                if 'ngal' in custom:
                    cl.ngal = cl.ngalaxies
                elif 'N_gal' in custom:
                    cl.N_gal = cl.ngalaxies
                if 'redshift' in custom:
                    cl.redshift = cl.z
                if 'iau_name' in custom:
                    cl.get_iau_name(prefix=prefix)
                if 'redshift_code' in custom:
                    cl.redshift_code = cl.z_type
                if 'theta80' in custom:
                    cl.theta80 = cl.r_x(frac=0.8)
                if 'concentration' in custom:
                    cl.concentration = cl.r_x(frac=0.8)/cl.r_x(frac=0.2)
                if 'scatter' in custom:
                    cl.scatter = cl.cmr(fit=True,colour=cl.basis)[1]

                if img_src:
                    cl.imgURL = cl.img_name()


        if custom:
            values = custom
        if img_src:
            values.append('imgURL')

        if len(extras) > 0:
            print 'Adding %s at your risk!'%(str(extras))
            values = values + list(extras)
            

        output = output.split('.fit')[0]+'.fit'
        out = cat(self,output=output,explicit=values)
        if thesis:
            del cl.RA
            del cl.DEC
            del cl.cluster_z
            del cl.cz_type
            del cl.N_gal
            del cl.bgc
            del cl.scatter
            del cl.theta80
            del cl.concentration

        else:
            from clusters import cluster as _cl
            dummy_cl = _cl([],-1,dummy=True)
            for cl in self.clusters:
                cl.redshift = cl.r_method
                for var in custom:
                    if not var in dummy_cl.__dict__.keys() and not var in preserve:
                        try:
                            del cl.__dict__[var]
                        except:
                            pass

        print 'Top-level written to %s \n\n'%(output)

            
        if texify:
            from tools import fits_to_tex
            cf = {}

            if thesis:
                cf = {'RA':"3.5f",'DEC':"3.5f",'cluster_z':"1.3f",'N_gal':"5.0f",'bgc':"5.0f",'scatter':"1.3f",'theta80':"1.4f",'concentration':"1.3f"}
                labels = {'iau_name':'Name','N_gal':"""N_{\\rm gal}""",'theta80':"""\\theta_{80}""",'concentration':"""\\theta_{80}/\\theta_{20}"""}

                if 'agc1i' in values:
#                    labels['agc1i'] = """A_{gc[1,\rm i]}"""
                    cf['agc1i'] = "1.3f"
                if 'stellarmass' in values:
#                    labels['stellarmass'] = """\rm{M_{\star}^{cl}}}"""
                    cf['stellarmass'] = "1.2f"


                if 1:
                    cf['bgc'] = "1.3f"

            else:
                cf = {'RA':"3.5f",'DEC':"3.5f",'cluster_z':"1.3f",'N_gal':"5.0f",'bgc':"5.0f",'scatter':"1.3f",'theta80':"1.4f",'concentration':"1.3f"}
                labels = {'iau_name':'Name','N_gal':"""N_{\\rm gal}""",'theta80':"""\\theta_{80}""",'concentration':"""\\theta_{80}/\\theta_{20}""",'redshift':'cluster_z'}
                labels['redshift_code'] = 'cz_type'




#                values = ['ID','ra','dec','ra_bcg','dec_bcg','ngal','redshift','redshift_code','theta80','concentration','stellarmass']
#                cf = {'RA':"3.5f",'DEC':"3.5f",'cluster_z':"1.3f",'N_gal':"5.0f",'bgc':"5.0f",'scatter':"1.3f",'theta80':"1.4f",'concentration':"1.3f",'stellarmass':"1.3f"}
#                labels = {'iau_name':'Name','N_gal':"""N_{\\rm gal}""",'theta80':"""\\theta_{80}""",'concentration':"""\\theta_{80}/\\theta_{20}""",'stellarmass':"""\rm{M_{\star}^{cl}}} """}




            print values
            out = fits_to_tex(output,order=values,custom_format=cf,labels=labels,longtable=True,landscape=True)
#            out = fits_to_tex(output,order=values,top=10,custom_format=cf)
#                cf = {'RA':"3.5f",'DEC':"3.5f",'bgc':"5.0f",'cluster_z':"1.3f",'b_gc':"d"}
#            out = fits_to_tex(output,order=values,top=10,custom_format=cf)



            return out

    def to_hdf5(self,output='out.hdf5',wrap=False,zip_level=6,updates=True,toplevel_fits=False,mags='auto',img_src=True):
        import h5py
        self.wrap_toggle(mode='observed')
#        mags = ['u','g','r','i','z']
        if mags == 'auto':
            mags = [self.selection['detection']['sequence'][0][0]]+[self.selection['colours'][c]['magnitude'] for c in self.selection['detection']['sequence']]


        status = None
        if updates:
            from dataIO import Progress
            print '\n\nConverting to HDF5'
            status = Progress(len(self.clusters))

        tree = h5py.File(output,'w')
        for cl in self.clusters:
            name = 'Cluster%d'%(indexOf(self.clusters,cl))
            cluster = tree.create_group(name)
            cx = degrees(cl.x)
            if wrap == False and cx < 0:
                cx = cx + 360.0
            cluster['ID'] = abs(cl.dna)
            cluster['ra'] = cx
#            cluster['ra'] = degrees(cx)
            cluster['dec'] = degrees(cl.y)

            bx = degrees(cl.bcg.x)
            if wrap == False and bx < 0:
                bx = bx + 360.0
#            cluster['ra_bcg'] = degrees(bx)
            cluster['ra_bcg'] = bx
            cluster['dec_bcg'] = degrees(cl.bcg.y)
            cluster['ngal'] = cl.ngalaxies
            cluster['redshift'] = cl.z
            cluster['redshift_code'] = cl.z_type

            cluster['theta80'] = cl.r_x(frac=0.8)
            cluster['concentration'] = cl.r_x(frac=0.8)/cl.r_x(frac=0.2)

            richness = tree[name].create_group('Richness')
            
            try:
                _agc = numpy.log10(cl.rich_data['bgc']['Agc'])
#                numpy.log10(cl.rich_data['bgcR_r80']['Agc'])
            except:
                _agc = numpy.nan
            if numpy.isnan(_agc) or numpy.isinf(_agc):
              rc = numpy.nan  
            elif _agc < -3.25:
                rc = 1
            elif (_agc >= -3.25) and (_agc < -2.625):
                rc = 2
            elif (_agc >= -2.625) and (_agc < -2.25):
                rc = 3
            elif _agc >= -2.25:
                rc = 4
            else:
                print 'Could not classify:'
                print _agc,self.dna
                raise SystemExit(0)
                
                               
            richness['Class'] = rc

            try:
                rd = cl.rich_data
            except:
                cl.rich_data = {}
                categories = ['bgc','bgcR','bgcR_r80','bgc_r80']
                values = ['Agc','N_b','N_t','bgc']
                for key in categories:
                    cl.rich_data[key] = {}
                    for sub_key in values:
                        cl.rich_data[key][sub_key] = numpy.nan
            

            types = {}
            sub_types = {}
            types['bgc'] = 'R500'
            types['bgcR'] = 'R500'
            types['bgcR_r80'] = 'theta80'
            types['bgc_r80'] = 'theta80'

            sub_types['bgc'] = 'MagLim'
            sub_types['bgcR'] = 'RedSequence'
            sub_types['bgcR_r80'] = 'RedSequence'
            sub_types['bgc_r80'] = 'MagLim'

            r500 = tree[name]['Richness'].create_group('R500')
            t80 = tree[name]['Richness'].create_group('theta80')

            r500.create_group('MagLim')
            r500.create_group('RedSequence')

            t80.create_group('MagLim')
            t80.create_group('RedSequence')



            for rclass in cl.rich_data.keys():
#                tree[name]['Richness'].create_group(rclass)
                for key in cl.rich_data[rclass].keys():
                    tree[name]['Richness'][types[rclass]][sub_types[rclass]][key] = cl.rich_data[rclass][key]

            galaxies = tree[name].create_group('Galaxies')
            for g in cl.galaxies:
                gname = 'Galaxy%d'%(indexOf(cl.galaxies,g))
                gal = tree[name]['Galaxies'].create_group(gname)
                gx = degrees(g.x)
                if wrap == False and gx < 0:
                    gx = gx + 360.0
                tree[name]['Galaxies'][gname]['ra'] = gx
                tree[name]['Galaxies'][gname]['dec'] = degrees(g.y)
                gal['objID'] = g.objID
                try:
                    gal['DR7_objID'] = g.extras['DR7objID']
                except:
                    pass
                for band in mags:
                    gal[band] = g.magnitudes[band]

                types = ['specz','2slaq','wiggleZ']
                try:
                    has_spec = numpy.array([not numpy.isnan(g.extras[s_type]) for s_type in types])
                except:
                    has_spec = [False]
                if True in has_spec:
                    gal['specz'] = g.extras[types[indexOf(has_spec,True)]]
                    gal['specz_source'] = types[indexOf(has_spec,True)]
                else:
                    gal['specz'] = g.specz
                    gal['specz_source'] = 'None'

                types = ['photoz','hyperz']
                try:
                    has_spec = numpy.array([not numpy.isnan(g.extras[s_type]) for s_type in types])
                except:
                    has_spec = [False]
                if True in has_spec:
                    gal['photoz'] = g.extras[types[indexOf(has_spec,True)]]
                    gal['photoz_source'] = types[indexOf(has_spec,True)]
                else:
                    gal['photoz'] = g.photoz
                    gal['photoz_source'] = 'None'
                
            if status:
                status.update(0)


                
        tree.close()
        print 'Detection written to %s\n\n'%(output)

        if toplevel_fits:
            from tools import cat
            for cl in self.clusters:
                cl.r_method = cl.redshift
                cl.ID = abs(cl.dna)
                if cl.x < 0:
                    cl.ra = degrees(cl.x) + 360
                else:
                    cl.ra = degrees(cl.x)
                cl.dec = degrees(cl.y)
                if cl.bcg.x < 0:
                    cl.ra_bcg = degrees(cl.bcg.x) + 360
                else:
                    cl.ra_bcg = degrees(cl.bcg.x)
                cl.dec_bcg = degrees(cl.bcg.y)
                cl.ngal = cl.ngalaxies
                cl.redshift = cl.z
                cl.redshift_code = cl.z_type
                cl.theta80 = cl.r_x(frac=0.8)
                cl.concentration = cl.r_x(frac=0.8)/cl.r_x(frac=0.2)
                cl.imgURL = cl.image
            values = ['ID','ra','dec','ra_bcg','dec_bcg','ngal','redshift','redshift_code','theta80','concentration']
            if img_src:
                values.append('imgURL')
            output = output.split('.hdf5')[0]+'.fit'
            out = cat(self,output=output,explicit=values)
            for cl in self.clusters:
                cl.redshift = cl.r_method
            print 'Top-level written to %s \n\n'%(output)
        self.wrap_toggle()

    def search(self,keyword='ngal'):
        from tools import search as sub_search
        results = []
        for key in self.selection.keys():
            if keyword in key:
                results.append('\%s\%s'%(key,self.selection[key]))
            if type(self.selection[key]) == type({}):
                path = '%s'%(key)
                sub_results = sub_search(self.selection[key],keyword,path=path)
                results = results + sub_results
        if len(results) > 0:
            print '\n'
            for result in results:
                print result
            print '\n'
        else:
            print 'No match'
        return None

    def sdssimg_parameters(self,server='http://cas.sdss.org/DR7/en/tools/search/x_sql.asp',batches=200,obj_override=None,verbose=False):
        if verbose: print 'Fetching run,rerun,camcol,field'
        from urllib import urlencode,urlopen
        import numpy
        context = 'PhotoObjAll'
        if '82' in self.selection['io']['dataset']:
            context = 'Stripe82..PhotoObjAll'
            server = 'http://cas.sdss.org/stripe82/en/tools/search/x_sql.asp'
        
        source = self.cluster_galaxies
        if obj_override != None:
            source = obj_override


        gid = numpy.array([g.objID for g in source])
        source = numpy.array(source)


        queries = []
        gals = iter(self.cluster_galaxies)
        query = 'select p.objID, p.run,p.rerun,p.camcol,p.field from %s as p where (p.objID = %d)'%(context,gals.next().objID)
        i = 1
        while 1:
            try:
                if i >= batches:
                    queries.append(query)
                    query = 'select p.objID, p.run,p.rerun,p.camcol,p.field from %s as p where (p.objID = %d)'%(context,gals.next().objID)
                    i = 0
                    i+=1
                else:
                    query = query + 'or (p.objID = %d) '%(gals.next().objID)
                    i+=1
            except:
                queries.append(query)
                break
        if verbose: print 'There are %d query blocks, %d objects'%(len(queries),len(source))
        results = []
        proc = []
        complete = []

        for q in queries:
            query_code = urlencode({'cmd':q,'format':'csv'})
            queryURL = server+'?'+query_code
            result = list(unique(urlopen(queryURL).read().splitlines()[1:]))
            results = results + result
            proc.append(result)
            complete.append(result)
            if len(results) > 55 and len(proc) > 55:
                print 'Close to submission limit, waiting for 60 seconds before commencing'
                import time
                time.sleep(60)
                print 'Commencing, currently on block %d/%d'%(len(complete),len(queries))
                proc = []


        if verbose: print 'Done generating, now sorting'
        for r in results:
            try:
                id,run,rerun,camcol,field = numpy.array(r.split(','),dtype='int')
            except:
                print 'Error:'
                print r
                return
            gg = source[numpy.where(gid == id)[0]]
            for g in gg:
                g.run,g.rerun,g.camcol,g.field = run,rerun,camcol,field

        source = list(source)
        if verbose: print 'Fin'
    

    def post_proc(self,rfrac=200.00,cleanup=True,clean_am=False):
        """\

        Scans through clusters, and allocates the brightest of all galaxies, inc
        associate members as the bcg

        """
        if len(self.clusters) == 0:
            return
        self.wrap_toggle(mode='detected')

        try:
            ucls = self.uber_clusters[0]
        except:
            self.retrieve_associates(uber=True,rfrac=rfrac)
        try:
            aa = self.associated_cluster_galaxies
        except:
            self.retrieve_associates(uber=True,rfrac=rfrac)
        try:
            aa = self.associated_cluster_galaxies[0]
        except:
            self.retrieve_associates(uber=True,rfrac=rfrac)

        from clusters import cluster
        det_clusters = []
        for cl in self.clusters:
#            if cl == self.clusters[6]:
#                print 'here'
            if len(cl.associated_members) > 0:
                #may have a brighter BCG....
                ucl = cl.uber_cluster
                if cl.bcg.dna != cl.uber_cluster.bcg.dna:
                    #does have a brighter BCG - need to generate a new cluster
                    cl.old_bcg = cl.bcg
                    lookup = {}
                    for g in cl.associated_members:
                        lookup[g.objID] = g
#                    cl.bcg = cl.uber_cluster.bcg
                    cl.bcg = lookup[cl.uber_cluster.bcg.objID]

                    galaxies = list(cl.galaxies) + [cl.bcg]
                    _cl = cluster(galaxies,0,real=True)
                    _cl.associate([cl])
                    _cl = cl.transplant(_cl,exclude=['galaxies','bcg','ngalaxies','associated_members','uber_cluster'])
                    _cl.associated_members = cl.associated_members
                    _cl.uber_cluster = cl.uber_cluster
                    for _g in _cl.associated_members:
                        _g.parent = _cl

                    _cl.basis = deepcopy(cl.basis)
                    _cl.selection = deepcopy(cl.selection)
                    _cl.nfriends = deepcopy(cl.nfriends)
                    _cl.truncation_history = deepcopy(cl.truncation_history)
                    _cl.friends = deepcopy(cl.friends)
                    _cl.rich_data = deepcopy(cl.rich_data)
                    try:
                        _cl.N_t = deepcopy(cl.N_t)
                        _cl.N_b = deepcopy(cl.N_b)
                        _cl.Agc = deepcopy(cl.Agc)
                        _cl.bgc = deepcopy(cl.bgc)
                    except:
                        pass

                    del cl
                    det_clusters.append(_cl)
                else:
                    det_clusters.append(cl)
            else:
                det_clusters.append(cl)
                
        self.addclusters(det_clusters,verbose=False)

        if cleanup:
            self.uber_clean(clean_am=clean_am)
        return
        if cleanup:
            del self.uber_cluster_galaxies
            del self.uber_clusters
            if clean_am:
                del self.associated_cluster_galaxies
            
            for cl in self.clusters:
                try:
                    if clean_am:
                        del cl.associated_members
                    del cl.uber_cluster
                except:
                    pass
                    
                    

    def uber_clean(self,clean_am=True):
        """\
        removes all uber-cluster and associate member data from detection class
        (will help sizes when saving det to disk)
        """
        try:
            del self.uber_cluster_galaxies
        except:
            pass
        try:
            del self.uber_clusters
        except:
            pass
        if clean_am:
            del self.associated_cluster_galaxies
            
        for cl in self.clusters:
            try:
                if clean_am:
                    del cl.associated_members
                del cl.uber_cluster
            except:
                pass
        



    def uber_cluster_generator(self,rfrac=0.5,pass_det=False,clobber=True):
        from clusters import galaxies_within_r,cluster
        from dataIO import Progress
        self.wrap_toggle(mode='detected')
        try:
            uc = self.uber_clusters
            if clobber == False:
                print 'Uber clusters already calculated'
                return
        except:
            pass

        cl = self.clusters[0]
        try:
            members = cl.associated_members
        except:
            print 'Retrieving associated memberships'
            self.retrieve_associates(pp=True)
        
        uc = []
        uc_gals = []
        print 'Generating uber clusters'
        status = Progress(len(self.clusters))
        for c in self.clusters:
            uber = c.generate_uber_cluster(rfrac=rfrac)
            uc.append(uber)
            for g in uber.galaxies:
                if (g in c.galaxies) == False:
                    uc_gals.append(g)
            status.update(0)
        self.uber_clusters = uc
        self.uber_cluster_galaxies = uc_gals
        if pass_det:
            from detections import detections as det
            dummy = det(None,deepcopy(self.selection),dummy=True)
            dummy.addclusters(self.uber_clusters)
            return dummy

    def filter_from_external(self,file,keyword='id',keyword2='',verbose=True,return_id=False,do_filter=True,mode='clusters',exclude=False,uber=False):
        """\
        example usage:
        det = det.filter_from_external('../boss/targetted_clusters.fits',do_filter=True,return_id=False,exclude=False)
        ids = det.filter_from_external('../boss/targetted_clusters.fits',do_filter=False,return_id=True)
        filtered = det2.filter_from_external('priority_supplemental.fit',keyword='gal_objID',exclude=False,mode='galaxy',uber=True,keyword2='clusterID')

        """
        import numpy
        from detections import detections as det
        dummy = det(None,deepcopy(self.selection),dummy=True)
        id2 = None
        if '.fit' in file:
            import fits_table as fits
            tbl = fits.fits_table(file)
            keys = tbl.names
            chosen = ''
            chosen2 = None
            for key in keys:
                options = [key,keyword.upper(),keyword.lower()]
                select = [keyword in key, keyword.upper() in key,keyword.lower() in key]
                if True in select:
                    chosen = options[indexOf(select,True)]
                    chosen = key
                if keyword2 != '':
                    options = [key,keyword2.upper(),keyword2.lower()]
                    select = [keyword2 in key, keyword2.upper() in key,keyword2.lower() in key]
                    if True in select:
                        chosen2 = options[indexOf(select,True)]
                        chosen2 = key



            id = tbl.get_column(chosen)
            if chosen2:
                id2 = tbl.get_column(chosen2)
        else:
            import numpy

        print 'id data read'
        if do_filter:
            if mode == 'clusters':
                cl_dna = numpy.array([c.dna for c in self.clusters])
                if exclude:
                    filter = numpy.array([(dna in id) == False for dna in cl_dna])
                else:
                    filter = numpy.array([dna in id for dna in cl_dna])
                self.clusters = numpy.array(self.clusters)
                clusters = self.clusters[filter]
                self.clusters = list(self.clusters)
                dummy.addclusters(clusters,verbose=False)
                return dummy
            elif 'gal' in mode:
                self.clusters = numpy.array(self.clusters)
                galaxies = self.cluster_galaxies
                if uber:
                    galaxies = galaxies + self.uber_cluster_galaxies
                objID = numpy.array([c.objID for c in galaxies])
                selected = []

                if id2 != None:
                    clfilter = numpy.array([cl.dna in id2 for cl in self.clusters])
                    cl_selected = self.clusters[clfilter]
                    for i in xrange(len(cl_selected)):
                        if uber:                            
                            _g = indexOf([g.objID for g in cl_selected[i].uber_cluster.galaxies],id[i])
                        else:
                            _g = indexOf([g.objID for g in cl_selected[i].galaxies],id[i])
                        selected.append(_g)
                else:
                    selected = numpy.array([g.objID in id for g in galaxies])
                    selected = galaxies[selected]
                    
                filter = numpy.array([g in selected for g in galaxies])

                if exclude:
                    filter = numpy.array([g not in selected for g in galaxies])

                galaxies = numpy.array(galaxies)
                galaxies = galaxies[filter]
                self.cluster_galaxies = list(self.cluster_galaxies)
                return galaxies
                


        elif return_id:
            return id


    def hectospec_measurements(self,plot=False):
        import tools
        md_num = int(self.selection['io']['dataset'].split('md0')[-1][0])
        return tools.hectospec_measurements(md_num,det=self,plot=plot,survey='mds')
    


    def hectospec_targetlist(self,ra=None,dec=None,outname='auto',genimg=True,rthresh=2.0,sloan_spec=None):
        import os
        import tools
        md_num = int(self.selection['io']['dataset'].split('md0')[-1][0])
        base_url = 'http://astro.dur.ac.uk/~dmurphy/orca/data/md0%d_hs/'%(md_num)
        measured = self.hectospec_measurements(plot=False)

        if len(measured) > 0:
            self.compare(measured,name='measured',type='gal')
        if sloan_spec and os.path.isfile(sloan_spec):
            print 'Reading in SDSS spec data - this may take some time...'
            self.compare(sloan_spec,name='sloan_spec',type='gal')
            print '...done'
        if ra and ':' in ra:
            ra = tools.ra_to_deg(ra)
        
        if dec and ':' in dec:
            dec = tools.dec_to_deg(dec)

        if ra and dec:
            d2 = det.find(ra,dec,r=60)        ##returns a list of clusters matched within 1 degree of the above position
        else:
            d2 = self

        img = []
        p1 = []                    ## BCG (not strictly in r)
        p2 = []                    ## 2nd brightest member (not strictly in r)
        p3 = []                    ## non-dominant members
        p4 = []                    ## if outside target area
        rlim = 21.5
        blacklist = [cl for cl in d2.clusters]
        no_fibres = False
        if (self.cluster_galaxies[0].extras.has_key('r_fiber') == False):
            print 'Warning, no fiber-mags for this detection class:\nfudging with r-band data (this could be wrong!)'
            gals = self.cluster_galaxies
            try:
                gals += self.associated_cluster_galaxies
            except:
                pass

            for g in gals:
                g.extras['r_fiber'] = g.magnitudes['r']
            no_fibres = True

        already_measured = []
        if self.external_data.has_key('measured') or self.external_data.has_key('sloan_spec'):
            from clusters import galaxies_within_r
            ext_gals = []
            if self.external_data.has_key('measured'):
                ext_gals += list(self.external_data['measured']['data'])
            if self.external_data.has_key('sloan_spec'):
                ext_gals += list(self.external_data['sloan_spec']['data'])

            allgals = list(self.cluster_galaxies)
            try:
                allgals += self.associated_cluster_galaxies
            except:
                pass
            for g in allgals:
                if len(galaxies_within_r(g,ext_gals,r=radians(rthresh/(60*60.0)))) > 0:
                    already_measured.append(g)
            if len(already_measured) > 0:
                print 'Blacklisted %d already-measured targets from selection'%(len(already_measured))

        for cl in d2.clusters:
            projects = iter(['bcg','b2cg','bncg','acm','iz_zy'])
            project = projects.next()
            g = cl.bcg
            if g.extras['r_fiber'] < rlim:
                if (g in already_measured) == False:
                    p1.append(g)
                    if genimg:
                        img = mds_quickimg(cl,scale=0.5,project=project)
                        img = img.split('/')[-1];g.url = base_url+project+'/%s'%(img)
            g = cl.nth_brightest(n=2)
            project = projects.next()
            if g.extras['r_fiber'] < rlim:
                if (g in already_measured) == False:
                    p2.append(g)
                    if genimg:
                        img = mds_quickimg(cl,scale=0.5,project=project,gal=g)
                        img = img.split('/')[-1];g.url = base_url+project+'/%s'%(img)
            jj = arange(3,len(cl.galaxies)+1)
            project = projects.next()
            for j in jj:
                g = cl.nth_brightest(n=j)
                if g.extras['r_fiber'] < rlim:
                    if (g in already_measured) == False:
                        p3.append(g)
                        if genimg:
                            img = mds_quickimg(cl,scale=0.5,project=project,gal=g)
                            img = img.split('/')[-1];g.url = base_url+project+'/%s'%(img)

            project = projects.next()
            if len(cl.associates['dna']) > 0:
                try:
                    asm = cl.associated_members
                except:
                    continue
                for g in cl.associated_members:
                    if g.extras['r_fiber'] < rlim and (cl in self.cluster_galaxies == False):
                        if (g in already_measured) == False:
                            p3.append(g)
                            if genimg:
                                img = mds_quickimg(cl,scale=0.5,project=project,gal=g)
                                img = img.split('/')[-1];g.url = base_url+project+'/%s'%(img)
                    
        project = 'iz_zy'
        basis = numpy.array([cl.basis for cl in self.clusters])
        filter = numpy.array([cl.basis == 'iz' and cl.dna not in blacklist for cl in self.clusters])
        self.clusters = numpy.array(self.clusters)
        izcl = self.clusters[filter]

        izcl = []

        for cl in izcl:
            g = cl.bcg
            if g.extras['r_fiber'] < rlim:
                if (g in already_measured) == False:
                    p4.append(g)
                    if genimg:
                        img = mds_quickimg(cl,scale=0.5,project=project)
                        img = img.split('/')[-1];g.url = base_url+project+'/%s'%(img)

        junk = """\
        project = 'full_field'
        if md_num == 8:
            for m in self.external_data['measured']['data']:
                blacklist += list(clusters_within_r(m,self.clusters,r=radians(10.0/60.0)))

        filter = numpy.array([cl.basis == 'ri' and cl not in blacklist for cl in self.clusters])
        remainder = self.clusters[filter]
        for cl in remainder:
            g = cl.bcg
            if g.extras['r_fiber'] < rlim:
                p4.append(g)
                if genimg:
                    img = mds_quickimg(cl,scale=0.5,project=project)
                    img = img.split('/')[-1];g.url = base_url+project+'/%s'%(img)

        """
        priorities = [p1,p2,p3,p4]
        all_targets = []
        print 'There are :'
        for i in xrange(len(priorities)):
            print '- %d targets in priority %d list'%(len(priorities[i]),i+1)
            for cl in priorities[i]:
                cl.priority = i+1
            all_targets = all_targets + list(priorities[i])

        data = {}
        data['ra'] = [degrees(g.x) for g in all_targets]
        data['dec'] = [degrees(g.y) for g in all_targets]
        data['r_fiber'] = [g.extras['r_fiber'] for g in all_targets]
        data['clusterID'] = [g.parent.dna for g in all_targets]
        data['objID'] = [g.objID for g in all_targets]
        data['priority'] = [g.priority for g in all_targets]
        order = ['objID','ra','dec','r_fiber','priority','clusterID']
        if genimg:
            data['imgURL'] = [g.url for g in all_targets]
            order.append('imgURL')


        if outname == 'auto':
            outname = 'md0%d_hs.fits'%(md_num)
        
        from tools import dict_to_fits
        dict_to_fits(data,outname,order=order)
        if no_fibres:
            print 'Now cleaning up r_fiber from extras'
            for g in gals:
                if g.extras.has_key('r_fiber'):
                    del g.extras['r_fiber']
            print 'Done'
        



    def find(self,ra=None,dec=None,r=1,units='arcmin',plot=False,buffer=1.5,verbose=False):
        """\
        det.find(self,ra=None,dec=None,r=1,units='arcmin',plot=False,buffer=1.5,verbose=False) => detection

        \t Search for clusters in this detection class, return new class containing recovered clusters
        
        ra         :  the right ascension to search at
        dec        :  the declination to search at
        r          :  the distance to look out to...
        units      :  ...in these units
        plot       :  display the results
        buffer     :  add a buffer to the limits delineating the search radius
        verbose    :  write status to screen
        """

        from tools import ra_to_deg,dec_to_deg
        from clusters import clusters_within_r

        r_input = deepcopy(r)

        if units in ['d','deg']:
            units = 'degrees'
        elif units in ['am','arcmin','min']:
            units = 'arcmins'
        elif units in ['as','arcsec','sec']:
            units = 'arcseconds'
        elif units in ['rad','rads']:
            units = 'radians'

        convert = {'degrees':1.0,'arcmins':60.0,'arcseconds':3600.0,'radians':None}
        if not units in convert.keys():
            print 'Unit %s not recognised'%(units)
            return None
        
        if not ra and not dec:
            return None
        if type(ra) == type(''):
            delimiters = [':',' ']
            delimiter = None
            for d in delimiters:
                if d in ra:
                    delimiter = d
            if not delimiter and ('h' in ra) and ('m' in ra):
                ra = ra.replace('h',':')
                ra = ra.replace('m',':')
                ra = ra.replace('s',':')
                delimiter = ':'
                print '...'
            if not delimiter:
                return None
            _ra = radians(ra_to_deg(ra,delim=delimiter))
        else:
            #assume degrees..
            _ra = radians(ra)

        if type(dec) == type(''):
            delimiters = [':',' ']
            delimiter = None
            for d in delimiters:
                if d in dec:
                    delimiter = d
            if not delimiter and ('d' in dec) and ('m' in dec):
                dec = dec.replace('h',':')
                dec = dec.replace('m',':')
                dec = dec.replace('s',':')
                delimiter = ':'
                print '...'
            if not delimiter:
                return None
            _dec = radians(dec_to_deg(dec,delim=delimiter))
        else:
            #assume degrees..
            _dec = radians(dec)
            
        
        if units != 'radians':
            _r = radians(r/convert[units])
#            print 'Converted %d %s to %f radians'%(r,units,_r)
            r = _r

        #
        self.wrap_toggle(mode='observed')
        #

#        print _ra,_dec
        nearby = clusters_within_r([_ra,_dec],self.clusters,r=r)
        if len(nearby) > 0:
            output = self.clone()
            from clusters import galaxy,cluster
            centre = galaxy(0,_ra,_dec,dummy=True)
            centre = cluster([centre],0,real=False,dummy=True)
            centre.extent = r
            centre.parent = output
            output.external_data = {}
            output.external_data['location'] = {'data':[centre], 'colour':'yellow', 'shape':'circle', 'points':False, 'y':[0.0], 'x':[0.0], 'type':'cluster'}
            output.addclusters(nearby,verbose=False)
            output.selection['limits']['xmin'] = max(0.0,_ra - buffer*r)
            output.selection['limits']['xmax'] = min(_ra + buffer*r,radians(360.0))
            output.selection['limits']['ymin'] = _dec - buffer*r
            output.selection['limits']['ymax'] = _dec + buffer*r
            if plot:
                output.plot(interact=True,wrap_mode='observed')

            if verbose:
                print '@ %s, %s found %d ORCA clusters within %s %s:'%(ra,dec,len(output.clusters),str(r_input),units)
                for cl in output.clusters:
                    print '  -- id=%d @ %3.4f, %3.4f'%(cl.dna,degrees(cl.x),degrees(cl.y))
                print ' '
            return output
        


        #
        #self.wrap_restore()
        #                         

        return []

            


    def save(self,filename):
        from dataIO import save
        save(self,filename)

    def clone(self,full=False):
        from copy import deepcopy
        from detections import detections
        if full == True:
            _det = deepcopy(self)
        else:
            _det = detections(None,deepcopy(self.selection),dummy=True)
        return _det

    def colour_dna(self):
        colours = self.selection['detection']['sequence']
        string = ''
        sel = self.selection['colours']
        for c in colours:
            s = '%1.8f %1.8f %1.8f'%(sel[c]['intercept'],sel[c]['slope'],sel[c]['width'])
            string = string + s
        dna = abs(hash(string))
        self.selection_dna = dna
        return dna


    def pp_richness(self,src_data=None,job_server=None,ncpus=16,uber=False,cleanup=True,clean_am=False,rfrac=0.8,cosma=False):
        import numpy
        gen_sd = False
        if src_data:
            self.src_data = src_data
            gen_sd = False
        try:
            sd = self.src_data
            if not sd:
                gen_sd = True
        except AttributeError:
            print 'Making src data from scratch'
            gen_sd = True

        if gen_sd:
            cdna = {}
            cdna_sel = {}
            self.src_data = {}
            for cl in self.clusters:
                cl_dna = cl.colour_dna()
                if cdna.has_key(cl_dna):
                    cdna[cl_dna].append(cl)
                else:
                    cdna[cl_dna] = [cl]
                cdna_sel[cl_dna] = cl.selection['colours']
            dets = []
            total_det = self.clone()
            for c in total_det.selection['detection']['sequence']:
                total_det.selection['colours'][c]['width'] = 100.0
            total_det.get_input()
            self.src_data['total'] = total_det.cuts


            colour_dna = cdna.keys()
            for id in colour_dna:
                _det = self.clone()
                for c in _det.selection['detection']['sequence']:
                    _det.selection['colours'][c]['width'] = cdna_sel[id][c]['width']
                    _det.selection['colours'][c]['intercept'] = cdna_sel[id][c]['intercept']
                    _det.selection['colours'][c]['slope'] = cdna_sel[id][c]['slope']
                dets.append(_det)
            args = [{'source[i]':None},{'perfect':False}]

            import detections as classes
            import detections
            import fits_table as fits
            from dataIO import save,load,pp_run,selection_fn,slice_data,make_mask,make_detections,getData,sel_edit,sel_edit2,Progress
            run_function = 'get_input'
            dependent_functions = (selection_fn,slice_data,make_mask,make_detections,classes,getData,fits,detections)
            modules = ('clusters','numpy','parallel_qhull','operator','fits_table','operator','detections')
            dets[0].get_input()
            detections,job_server = pp_run(list(dets),run_function,args,dependent_functions,modules,ncpus=ncpus,idle_delay=360.0,server=job_server,verbose=False)

            print 'detections made'
            for d,id in zip(detections,cdna):
                self.src_data[id] = d.cuts


        clusters = numpy.array(self.clusters)
        cdna = numpy.array([_cl.colour_dna() for _cl in self.clusters])
        ucdna = numpy.array(unique(cdna))

####################
#        for dna in ucdna:
#            self.src_data[dna] = [dna for i in arange(5)]
####################

        for dna in cdna:
            if dna not in self.src_data.keys():
                print 'Missing colour_dna record %d'%(dna)

        nbar = max(int(float(len(cdna))/float(ncpus)),ncpus)
        print 'Aiming for %d clusters / detection class'%(nbar)
        dna_ref = {}
        for dna in ucdna:
            dna_ref[dna] = []
        for cl in clusters:
            dna_ref[cl.colour_dna()].append(cl)

        keys = numpy.array(dna_ref.keys())
        counts = numpy.array([len(dna_ref[k]) for k in keys])
        counts_srt = list(deepcopy(counts))
        counts_srt.sort()
        counts_srt.reverse()

        if uber:
            try:
                ucls = self.uber_clusters[0]
            except:
                self.retrieve_associates(uber=True,rfrac=rfrac)
            try:
                aa = self.associated_cluster_galaxies
            except:
                self.retrieve_associates(uber=True,rfrac=rfrac)
            try:
                aa = self.associated_cluster_galaxies[0]
            except:
                self.retrieve_associates(uber=True,rfrac=rfrac)
            
        dets = []
        in_keys = []
        while len(in_keys) != len(ucdna):
            det = self.clone()
            clusters = []
            for c in counts_srt:
                _keys = keys[numpy.where(counts == c)[0]]
                for key in _keys:
                    if not key in in_keys:
                        clusters = clusters + dna_ref[key]
                        in_keys.append(key)
                        if len(clusters) > nbar or len(in_keys) == len(ucdna):
                            det.addclusters(clusters,verbose=False)
                            det_dna = unique(numpy.array([_cl.colour_dna() for _cl in det.clusters]))
                            det.src_data = {'total':self.src_data['total']}
                            for dna in det_dna:
                                det.src_data[dna] = self.src_data[dna]
                            dets.append(det)
                            det = self.clone()
                            clusters = []
                        if len(in_keys) == len(ucdna):
                            break

        #sanity check:
        all_dna = []
        for d in dets:
            for dna in d.src_data.keys():
                if not (type(dna) == type('')):
                    all_dna.append(dna)

        if len(all_dna) == len(ucdna) and len(dets) > 0:
            print 'Passed sanity check, have %d detection classes'%(len(dets))
            if len(dets) > 1:
                rich_dets = pp_richness_launcher(dets,self,job_server=job_server,construct=False,uber=uber,rfrac=rfrac)
            else:
                rich_dets = [cluster_simple_richness(dets[0],do_all=True,uber=uber,rfrac=rfrac)]
            all_clusters = []
            for rd in rich_dets:
                all_clusters = all_clusters + list(rd.clusters)

            self.addclusters(all_clusters,verbose=False)
            for cl in self.clusters:
                cl.rich_data = deepcopy(cl.rich_data)
#            del self.src_data
            print 'Cluster richnesses calculated'
        if uber and cleanup:
            self.uber_clean(clean_am=clean_am)


    def cluster_simple_richness(self,pp=False,status_update=True,job_server=None,uber=False,rfrac=0.8,bgcR=False,r80=False,bcg_centre=False,do_all=False):
        from dataIO import Progress,pp_run
        selection_backup = deepcopy(self.selection)
        clusters = self.clusters
            
        if uber:
            try:
                clusters = self.uber_clusters
                cl = clusters[0]
            except:
                self.uber_cluster_generator(rfrac=rfrac)
                clusters = self.uber_clusters
        if status_update:
            status = Progress(len(self.clusters))
        for c in self.clusters:
            _c = c
            if uber:
                _c = c.uber_cluster
                
            _c.bgc_simple(bgcR=bgcR,r80=r80,bcg_centre=bcg_centre,do_all=do_all)
            c.rich_data = deepcopy(_c.rich_data)
            
            try:
                if c.rich_data['bgc']['Agc'] >= 0:
                    _agc = numpy.log10(c.rich_data['bgc']['Agc'])
                else:
                    _agc = numpy.nan
            except:
                _agc = numpy.nan
            if numpy.isnan(_agc) or numpy.isinf(_agc):
              rc = numpy.nan  
            elif _agc < -3.25:
                rc = 1
            elif (_agc >= -3.25) and (_agc < -2.625):
                rc = 2
            elif (_agc >= -2.625) and (_agc < -2.25):
                rc = 3
            elif _agc >= -2.25:
                rc = 4
            else:
                print 'Could not classify:'
                print _agc,self.dna
            c.richness_class = rc
            

            if status_update:
                status.update(0)
        
        
        

    def cluster_richness(self,pp=False,status_update=True,job_server=None,uber=False,rfrac=0.5,bgcR=False,r80=False,bcg_centre=False,do_all=False):
        from dataIO import Progress,pp_run
        selection_backup = deepcopy(self.selection)
        clusters = self.clusters
        if uber:
            try:
                clusters = self.uber_clusters
                cl = clusters[0]
            except:
                self.generate_uber_clusters(rfrac=rfrac)
                clusters = self.uber_clusters
        if pp:
            cl_list = {}
            for c in clusters:
                try:
                    cl_list[c.sub_id].append(c)
                except KeyError:
                    cl_list[c.sub_id] = [c]
            from detections import detections
            dets = []
            for key in cl_list.keys():
                _det = detections(None,deepcopy(self.selection),dummy=True)
                _det.addclusters(cl_list[key],verbose=False)
                dets.append(_det)
            args = [{'source[i]':None},{'bcgR':bgcR},{'r80':r80},{'bcg_centre':bcg_centre},{'do_all':do_all}]
            ncpus = 16
            dependent_functions = ()
            modules = ('science','parallel_detect','clusters','numpy','parallel_qhull','parallel_vertices','vperc','time','detections')
            run_function = 'cluster_richness'
            completed,job_server = pp_run(list(dets),run_function,args,dependent_functions,modules,ncpus=ncpus,idle_delay=360.0,verbose=True,server=job_server)
            clusters = []
            for d in completed:
                clusters = clusters + list(d.clusters)
            if uber:
                self.uber_clusters = clusters
                for c in self.uber_clusters:
                    c.parent = self
            else:
                self.addclusters(clusters,verbose=False)
            return

        else:
            if status_update: status = Progress(len(self.clusters))
            if bgcR or do_all:
                from tools import which_limits,which_subindex

                sids = []
                for c in clusters:
                    _s = which_subindex(c.x,c.y,wrap=True)
                    if type(_s) == type([]):
                        sids = sids + list(_s)
                    else:
                        sids.append(_s)
                sids = list(unique(sids))


                if bgcR or do_all:               
                    self.src_data = {}
                    for key in self.selection['colours'].keys():
                        self.selection['colours'][key]['width'] = 100.0
                
                    for s in sids:
                        self.selection['detection']['subset_index'] = s
                        self.selection['limits']['xmin'],self.selection['limits']['xmax'] = which_limits(s)
                        self.cuts = None
                        self.source_data = None
                        self.get_input()
                        self.src_data[s] = self.cuts
                    
                    print 'bgcR prep done'
            for c in self.clusters:
                c.bgc_richness(bgcR=bgcR,r80=r80,bcg_centre=bcg_centre,do_all=do_all)
                if status_update:status.update(0)
            try:
                del self.subset_cuts
            except:
                pass
        self.selection = selection_backup

        return

    def cz_scaling(self,cmap='auto'):
        import matplotlib as mpl
        for cl in self.clusters:
            cl.redshift()
        c = ['#0000FF','#FF0000']

        if type(cmap) == type('str') and cmap != 'auto':
            cmaps = {}
            cmaps['binary'] = plt.cm.binary
            cmaps['autumn'] = plt.cm.autumn

            try:
                cmap = cmaps[cmap]
            except KeyError:
                print 'cmap %s not recognise, reverting to auto'%(cmap)
                cmap = 'auto'
        elif type(cmap) != type('str'):
            if not cmap in plt.cm.__dict__.values():
                print 'Supplied cmap is (possibly?) an object, but not a recognised plt.cmap. Reverting to auto...'
                cmap = 'auto'
        if cmap == 'auto':
            cmap=mpl.colors.LinearSegmentedColormap.from_list('mycm',c,N=512)
                

        cz = [cl.z for cl in self.clusters]
        zmax = max(cz)
        zmin = min(cz)
        scaling = []
        for cl in self.clusters:
            scale = (cl.z-zmin)/(zmax-zmin)
#            scale = 1.0-scale
            cl.cz_scale = cmap(scale)
        
            


    def learn(self,minbin=[1,1,1,1],plot=False,mode='all',s=[1,1,1,1],weight=False,plot_bins=False,bin_samp_freq=[0.05,0.1,0.1,0.1],colours=['gr','ri'],der=0,cross_basis=False,sel_update=False):
        """\
        det.learn()

        \t Use the cluster data to learn about the evolution of sequence observables.
        
        minbin     :  the minimum # of clusters in a bin for it to be considered: [5,5,5] is safe
        plot       :  plot the results? (g-r = blue, r-i = red, i-z = green)
        mode       :  what are we correlating? [cm20-slope, slope-z, scatter-z]
        pp         :  run checks in parallel
        ncpu       :  #cpu to run in parallel

        """


        from scipy import interpolate
        from tools import bin_hist
        try:
            basis = numpy.array([cl.basis for cl in self.clusters])
        except:
            print 'Clusters have no basis data - setting cross_basis = True!'
            basis = numpy.ones(len(self.clusters))
            cross_basis = True
        self.clusters = numpy.array(self.clusters)

        smooth_auto = False
        if s[0] == 'a' or s[0] == 'auto':
            print 'automatic smoothing condition'
            smooth_auto = True

        if mode == 'all':
            print 'WARNING - the <<all>> mode does not work at the moment!'
            return None

        colour_wheel = iter(['b','r','g','magenta','orange','cyan'])
        colour_keys = {}
        colour_keys['gr'] = colour_wheel.next()
        colour_keys['ri'] = colour_wheel.next()
        colour_keys['iz'] = colour_wheel.next()
        for c in colours:
            if not c in colour_keys.keys():
                colour_keys[c] = colour_wheel.next()
                

        junk = """\
        try:
            fits = self.learning
            if not f.has_key('gr'):
                fits['gr'] = {}
            if not f.has_key('ri'):
                fits['ri'] = {}
        except:
            fits = {}
            fits['gr'] = {}
            fits['ri'] = {}
        """
        fits = {}
        for colour in colours:
            fits[colour] = {}


        colour1,colour2 = colours[0],colours[1]


        source_data = {}
        modes = ['z','slope','scatter','cm20']
        modes = {}
        modes['z'] = '[cl.z for cl in cls]'
        modes['slope'] = '[cl.cmr(fit=True)[0] for cl in cls]'
        modes['scatter'] = 'numpy.array([cl.cmr(colour=cl.basis,new_fit=True,magfilter=3,simple2sig=False)[1] for cl in cls])'
        modes['cm20'] = 'numpy.array([cl.get_cm20(colour=cl.basis) for cl in cls])'

        if cross_basis:       #means all clusters considered, not just those with matching basis (helps for i-z slope under {ri,iz})
            modes['slope'] = '[cl.cmr(fit=True,colour=colour)[0] for cl in cls]'
            modes['scatter'] = 'numpy.array([cl.cmr(colour=colour,new_fit=True,magfilter=3,simple2sig=False)[1] for cl in cls])'
            modes['cm20'] = 'numpy.array([cl.get_cm20(colour=colour) for cl in cls])'



        for colour in colours:
            source_data[colour] = {}
            f = (basis == colour)
            cls = self.clusters[f]
            if cross_basis:
                cls = self.clusters

            for _mode in modes.keys():
                if _mode in mode:
                    source_data[colour][_mode] = eval(modes[_mode])
                    
                  
        #z-slope
        if mode == 'z-slope':
            if plot:
                fig = plt.figure()
            for colour in colours:
                _sf = bin_samp_freq[0]
                if type(_sf) == type({}):
                    try:
                        _sf = _sf[colour]
                    except:
                        raise SystemExit('Sampling frequency not found for colour %s'%(colour))
                xr = numpy.arange(0,max(source_data[colour]['z']),_sf)
                out = bin_hist(source_data[colour]['z'],source_data[colour]['slope'],xr,colour='b',minbin=minbin[0],ls='--',plot=False,err=weight)
                wt = None
                if weight:
                    _z,_s,_err = out
                    wt = numpy.array(_err)
                    wt = 1.0/wt
                else:
                    _z,_s = out
                if smooth_auto == True:
                    _m = len(_z)
                    minm = _m-((2*_m)**0.5)
                    maxm = _m=((2*_m)**0.5)
                    dm = maxm-minm
                    s[0] = minm+(0.5*dm)
                    print 'Setting s[0] = %f'%(s[0])

                spline_zs = interpolate.splrep(_z,_s,s=s[0],w=wt)
                fits[colour]['z-slope'] = spline_zs
                if plot:
                    if plot_bins:
                        xr = numpy.arange(0,max(source_data[colour]['z']),_sf)
                        out = bin_hist(source_data[colour]['z'],source_data[colour]['slope'],xr,colour=colour_keys[colour],minbin=minbin[0],ls='--',plot=True,err=weight)
                    zrange = arange(min(_z),max(_z),0.01)
                    plt.plot(source_data[colour]['z'],source_data[colour]['slope'],color=colour_keys[colour],marker=',',ls='None')
                    yspline = interpolate.splev(zrange,fits[colour]['z-slope'],der=der)
                    plt.plot(zrange,yspline,color=colour_keys[colour],ls='-',marker='None')
                    if colour == colours[-1]:
                        ylim = plt.ylim()
                        plt.ylim(ylim[1],ylim[0])
                        plt.xlabel('redshift')
                        plt.ylabel('sequence slope')
                    
#            plt.show()
        #cm20-slope
        if mode == 'cm20-slope':
            if plot:
                fig = plt.figure()
            for colour in colours:
                _sf = bin_samp_freq[1]
                if type(_sf) == type({}):
                    try:
                        _sf = _sf[colour]
                    except:
                        raise SystemExit('Sampling frequency not found for colour %s'%(colour))

                xr = arange(min(source_data[colour]['cm20']),max(source_data[colour]['cm20']),_sf)
                out = bin_hist(source_data[colour]['cm20'],source_data[colour]['slope'],xr,minbin=minbin[1],ls='--',plot=False,err=weight)
                wt = None
                if weight:
                    _cm20,_slope,_err = out
                    wt = numpy.array(_err)
                    wt = 1.0/wt
                else:
                    _cm20,_slope  = out


                if smooth_auto == True:
                    _m = len(_cm20)
                    minm = _m-((2*_m)**0.5)
                    maxm = _m=((2*_m)**0.5)
                    dm = maxm-minm
                    s[1] = minm+(0.5*dm)
                    print 'Setting s[1] = %f'%(s[1])


                spline_cs = interpolate.splrep(_cm20,_slope,s=s[1],w=wt)
                fits[colour]['cm20-slope'] = spline_cs
                if plot:
                    if plot_bins:
                        xr = arange(min(source_data[colour]['cm20']),max(source_data[colour]['cm20']),_sf)
                        out = bin_hist(source_data[colour]['cm20'],source_data[colour]['slope'],xr,colour=colour_keys[colour],minbin=minbin[1],ls='--',plot=True,err=weight)
                    crange = arange(min(_cm20),max(_cm20),0.005)
                    plt.plot(source_data[colour]['cm20'],source_data[colour]['slope'],marker=',',color=colour_keys[colour],ls='None')
                    yspline = interpolate.splev(crange,fits[colour]['cm20-slope'],der=der)
                    plt.plot(crange,yspline,color=colour_keys[colour],ls='-',marker='None')
                    if colour == colours[-1]:
                        ylim = plt.ylim()
                        plt.ylim(ylim[1],ylim[0])
                        plt.xlabel('sequence normalisation')
                        plt.ylabel('sequence slope')


#            plt.show()


        #cm20-z
        if mode == 'cm20-z':
            if plot:
                fig = plt.figure()
            for colour in colours:
                _sf = bin_samp_freq[2]
                if type(_sf) == type({}):
                    try:
                        _sf = _sf[colour]
                    except:
                        raise SystemExit('Sampling frequency not found for colour %s'%(colour))

                xr = arange(min(source_data[colour]['cm20']),max(source_data[colour]['cm20']),_sf)
                out = bin_hist(source_data[colour]['cm20'],source_data[colour]['z'],xr,colour='b',minbin=minbin[2],ls='--',plot=False,err=weight)
                wt = None
                if weight:
                    _cm20,_z,_err = out
                    wt = numpy.array(_err)
                    wt = 1.0/wt
                else:
                    _cm20,_z = out

                if smooth_auto == True:
                    _m = len(_cm20)
                    minm = _m-((2*_m)**0.5)
                    maxm = _m=((2*_m)**0.5)
                    dm = maxm-minm
                    s[2] = minm+(0.5*dm)
                    print 'Setting s[2] = %f'%(s[2])


                spline_cz = interpolate.splrep(_cm20,_z,s=s[2],w=wt)
                fits[colour]['cm20-z'] = spline_cz

                if plot:
                    if plot_bins:
                        xr = arange(min(source_data[colour]['cm20']),max(source_data[colour]['cm20']),_sf)
                        out = bin_hist(source_data[colour]['cm20'],source_data[colour]['z'],xr,colour=colour_keys[colour],minbin=minbin[2],ls='--',plot=True,err=weight)
                    crange = arange(min(_cm20),max(_cm20),0.005)
                    plt.plot(source_data[colour]['cm20'],source_data[colour]['z'],marker=',',color=colour_keys[colour],ls='None')
                    yspline = interpolate.splev(crange,fits[colour]['cm20-z'],der=der)
                    plt.plot(crange,yspline,color=colour_keys[colour],ls='-',marker='None')


                    if colour == colours[-1]:
                        plt.xlabel('sequence normalisation')
                        plt.ylabel('redshift')

        #cm20-scatter
        if mode == 'cm20-scatter':
            if plot:
                fig = plt.figure()
            for colour in colours:
                _sf = bin_samp_freq[3]
                if type(_sf) == type({}):
                    try:
                        _sf = _sf[colour]
                    except:
                        raise SystemExit('Sampling frequency not found for colour %s'%(colour))

                xr = arange(min(source_data[colour]['cm20']),max(source_data[colour]['cm20']),_sf)
                out = bin_hist(source_data[colour]['cm20'],source_data[colour]['scatter'],xr,minbin=minbin[3],ls='--',plot=False,err=weight)

                wt = None
                if weight:
                    _cm20,_scatter,_err = out
                    wt = numpy.array(_err)
                    wt = 1.0/wt
                else:
                    _cm20,_sscatter  = out

                if smooth_auto == True:
                    _m = len(_cm20)
                    minm = _m-((2*_m)**0.5)
                    maxm = _m=((2*_m)**0.5)
                    dm = maxm-minm
                    s[3] = minm+(0.5*dm)
                    print 'Setting s[3] = %f'%(s[3])

                spline_cs = interpolate.splrep(_cm20,_scatter,s=s[3],w=wt)
                fits[colour]['cm20-scatter'] = spline_cs


                if plot:
                   if plot_bins:
                       xr = arange(min(source_data[colour]['cm20']),max(source_data[colour]['cm20']),_sf)
                       out = bin_hist(source_data[colour]['cm20'],source_data[colour]['scatter'],xr,colour=colour_keys[colour],minbin=minbin[3],ls='--',plot=True,err=weight)

                   crange = arange(min(_cm20),max(_cm20),0.005)
                   plt.plot(source_data[colour]['cm20'],source_data[colour]['scatter'],marker=',',color=colour_keys[colour],ls='None')
                   yspline = interpolate.splev(crange,fits[colour]['cm20-scatter'],der=der)
                   plt.plot(crange,yspline,color=colour_keys[colour],ls='-',marker='None')
                   if colour == colours[-1]:
                       plt.xlabel('sequence normalisation')
                       plt.ylabel('sequence scatter') 





#        _fits = {}
#        _fits[colour1] = fits['gr']
#        _fits[colour2] = fits['ri']
#        fits = _fits
        self.learning = fits
        if sel_update == True:
            for colour in fits.keys():
                if mode != 'all':
                    self.selection['colours'][colour]['fitfuncs'][mode] = fits[colour][mode]
        return fits


    def width(self,colour='default',cm20='default',der=0):
        from scipy import interpolate
        if colour == 'default':
            colour = self.selection['defaults']['colour']
        if cm20 == 'default':
            cm20 = self.selection['colours'][colour]['intercept']
        w = self.selection['colours'][colour]['width']
        if self.selection['colours'][colour]['fitfuncs'].has_key('cm20-scatter'):
            w = interpolate.splev(cm20,selection['colours'][colour]['fitfuncs']['cm20-scatter'],der=der)
            
        return w

    def slope(self,colour='default',cm20='default',der=0):
        from scipy import interpolate
        if colour == 'default':
            colour = self.selection['defaults']['colour']
        if cm20 == 'default':
            cm20 = self.selection['colours'][colour]['intercept']
        s = self.selection['colours'][colour]['slope']
        if self.selection['colours'][colour]['fitfuncs'].has_key('cm20-slope'):
            s = interpolate.splev(cm20,selection['colours'][colour]['fitfuncs']['cm20-slope'],der=der)
            
        return s

    def bcg_outliers(self,plot=False,sigma_clip=False,spine=False,ext_data=None):
        from numpy import std, mean, median
        import tools
        from scipy import interpolate
        det = self
        cls = det.clusters
        gr = numpy.array([cl.bcg.colours['gr'] for cl in cls])
        ri = numpy.array([cl.bcg.colours['ri'] for cl in cls])
        iz = numpy.array([cl.bcg.colours['iz'] for cl in cls])
        zz = [cl.z for cl in cls]
        mb=10
        zrange = arange(0,0.7,0.02)
        grspl = tools.bin_hist(zz,gr,zrange,err=True,spline=True,minbin=mb)[-1]
        rispl = tools.bin_hist(zz,ri,zrange,err=True,spline=True,minbin=mb)[-1]
        izspl = tools.bin_hist(zz,iz,zrange,err=True,spline=True,minbin=mb)[-1]
        gr_offset = numpy.array(interpolate.splev(zz,grspl,der=0))
        ri_offset = numpy.array(interpolate.splev(zz,rispl,der=0))
        iz_offset = numpy.array(interpolate.splev(zz,izspl,der=0))
        gr_adj = gr-gr_offset
        ri_adj = ri-ri_offset
        iz_adj = iz-iz_offset

        gr_std = std(gr_adj)
        ri_std = std(ri_adj)
        iz_std = std(iz_adj)

        if sigma_clip != False:
            if type(sigma_clip) != type(3):
                sigma_clip = 3
                
            filter_gr = numpy.array([c < 0 and abs(c) > sigma_clip*gr_std for c in gr_adj])
            filter_ri = numpy.array([c < 0 and abs(c) > sigma_clip*ri_std for c in ri_adj])
            filter_iz = numpy.array([c < 0 and abs(c) > sigma_clip*iz_std for c in iz_adj])

            filter_gr = numpy.array([abs(c) > sigma_clip*gr_std for c in gr_adj])
            filter_ri = numpy.array([abs(c) > sigma_clip*ri_std for c in gr_adj])
            filter_iz = numpy.array([abs(c) > sigma_clip*iz_std for c in gr_adj])

            gr_adj = gr_adj[filter_gr]
            ri_adj = ri_adj[filter_ri]
            iz_adj = iz_adj[filter_iz]
            
        



        if plot:
            plt.plot(gr_adj,zz,'b,')
            plt.plot(ri_adj+2,zz,'r,')
            plt.plot(iz_adj+4,zz,'g,')

#            plt.plot(gr,zz,'b,')
#            plt.plot(ri,zz,'r,')
#            plt.plot(iz,zz,'g,')




    def widen_selection(self,colour='gr',width=0.1,origin='base',new_input=False,plot=False):
        """\
        widens the selection and moves cm20 to straddle at half-width
        width is from the base (ie, the lower limit of the selection filter)
        unless origin is overridden
        """
        w = self.selection['colours'][colour]['width']
        if plot:
            self.plot_selection(colour=colour)
        w_ = width/2.0
        if w_ < w:
            print 'New width is smaller!'
            return
        dw = w_-w
        grad = self.selection['colours'][colour]['slope']
        from math import atan2,cos
        angle = atan2(-1.0*grad,-1.0)     #is this robust to gradient sign?
        cm20 = abs(dw/cos(angle))
        if origin == 'base':
            self.selection['colours'][colour]['intercept'] = self.selection['colours'][colour]['intercept'] + cm20
        self.selection['colours'][colour]['width'] = w_
        if plot:
            xmin,xmax = plt.xlim()
            self.plot_selection(colour=colour,plot_colour='g',passback=True)
            plt.xlim(xmin,xmax)
        if new_input:
            print 'Generating new set of galaxies'
            self.get_input()


    def plot_normalisation(self,passback=False):
        if passback == False:
            fig = plt.figure()
        grads,cm20s = [],[]
        for c in self.clusters:
            gg = c.cmr(fit=True)
            grads.append(gg[0])
            cm20s.append(gg[-1])

        plt.plot(grads,cm20s,'k.',ls='None')
        plt.xlabel('gradient')
        plt.ylabel('cm20')
        return
        
            
    def re_posess(self):
        self.addclusters(self.clusters,perfect=self.perfect_clusters,verbose=False)

    def retrieve_associates(self,pp=False,keep_cuts=False,status=True,job_server=None,uber=False,rfrac=0.5,gcross=True):
        from dataIO import Progress,pp_run
        from clusters import galaxy
        from detections import detections as det
        from copy import deepcopy
        dummy = det(None,deepcopy(self.selection),dummy=True)

        if pp:
            cl_list = {}
            for c in self.clusters:
                try:
                    cl_list[c.sub_id].append(c)
                except KeyError:
                    cl_list[c.sub_id] = [c]
            from detections import detections
            dets = []
            for key in cl_list.keys():
                _det = det(None,deepcopy(self.selection),dummy=True)
                _det.addclusters(cl_list[key],verbose=False)
                dets.append(_det)
            args = [{'source[i]':None},{'keep_cuts':keep_cuts}]
            ncpus = 16
            dependent_functions = ()
            modules = ('science','parallel_detect','clusters','numpy','parallel_qhull','parallel_vertices','vperc','time','detections')
            run_function = 'retrieve_associates'
            completed,job_server = pp_run(list(dets),run_function,args,dependent_functions,modules,ncpus=ncpus,idle_delay=360.0,verbose=False,server=job_server)
            clusters = []
            amgs = []
            for d in completed:
                clusters = clusters + list(d.clusters)
                amgs = amgs + list(d.associated_cluster_galaxies)


            self.addclusters(clusters,verbose=False)
            self.associated_cluster_galaxies = amgs


            if uber:
                self.uber_cluster_generator(rfrac=rfrac)
                self.addclusters(clusters,verbose=False)

            
            if grcoss:
                self.gcross()
                


            return




        for colour in dummy.selection['detection']['sequence']:
            dummy.selection['colours'][colour]['width'] = 100.0

        sid = False
        if 1:
            sid = self.clusters[0].sub_id
            sids = [c.sub_id for c in self.clusters]
            if len(unique(sids)) == 1:
                if sid == None and self.selection['detection'].has_key('subset_index') and self.selection['detection']['subset_index']  == -1:
                    sid = -1
                    
                dummy.selection['detection']['subset_index'] = sid
        else:
            pass

        try:
            dummy.cuts = self.src_data['total']
            print 'Using detection source data'
        except:
            dummy.get_input()
        uid = dummy.cuts['extras']['id']
        colours = dummy.selection['detection']['sequence']
        magnitudes = []
        for c in colours:
            magnitudes.append(dummy.selection['colours'][c]['magnitude'])

        if status:
            from dataIO import Progress
            status = Progress(len(self.clusters))
        
        dnas = None
        if dummy.cuts['extras']['id'] == None:
            dnas = [abs(hash('%1.8f %1.8f'%(float(dummy.cuts['x'][_i]),float(dummy.cuts['y'][_i])))) for _i in xrange(len(dummy.cuts['y']))]
        for cluster in self.clusters:
            cluster.retrieve_associates(dummy=dummy,dna_data=dnas)
            if status:
                status.update(0)
        amgs = []
        for cluster in self.clusters:
            amgs = amgs + cluster.associated_members
        self.associated_cluster_galaxies = amgs
        
        if uber:
            self.uber_cluster_generator(rfrac=rfrac)




    def maxBCG_completeness(self,nstep=10,factor=1,iterations=10,err=True,tolerance=5.0,r_override=False,maxBCG=None,plot=False):
        #loop through external data
        #look within maxBCG cluster radius 
        #are there >5 cluster galaxies residing in that radius?
        from clusters import galaxies_within_r,clusters_within_r
        from numpy import radians,arange,random,degrees,std,mean,median
        from clusters import cluster,galaxy
        from operator import countOf,indexOf
        from dataIO import Progress

        def counter(data,maxBCG,radii):
            for g in data:
                g.alloc = False
            data = numpy.array(data)
            counts = []

            for r in radii:
                count = 0
                filter = numpy.array([g.alloc == False for g in data])
                if len(data) > 0:
                    data = data[filter]
                else:
                    data = []
                for bcg in maxBCG:
                    nearby = numpy.array(galaxies_within_r(bcg,data,r=r))
                    if len(nearby) > 0:
                        dist = numpy.array([g.distance_from(xc=bcg.x,yc=bcg.y) for g in nearby])
                        select = (dist > r-dr)
                        nearby = nearby[select]
                        dist = dist[select]
                        filter = numpy.array([g.alloc == False for g in nearby])
                        nearby = nearby[filter]
                        dist = dist[select]
                        count = count + len(nearby)
                        for g in nearby:
                            g.alloc = True
                            g.sep = dist[indexOf(nearby,g)]
                counts.append(count)
            return counts


        #set up maxBCG
        if maxBCG == None:
            maxBCG = self.convert_from_external('../maxbcg.cat',wrap=True)
        self.maxBCG = maxBCG
        for c in self.maxBCG:
            c.parent = self
        rmax = max([c.extent for c in maxBCG])
        dr = rmax/float(nstep)
        radii = arange(dr,rmax+(dr/2.0),dr)

        #STAGE 1: BCG direct match - within <tolerance> arcsec

        detected_bcg = numpy.array([c.bcg for c in self.clusters])
        for c in self.clusters:
            c.bcg_match = False
            c.match_to = False
        for c in maxBCG:
            c.matches = []
            c.cg_bcg_match = False
            c.bcg_match = False
            c.stat_sig = False
            c.null_sig = False

        maxBCG = numpy.array(maxBCG)
        maxBCG_backup = maxBCG
        

        bcg_match = 0

        for bcg in maxBCG:
            filter = numpy.array([degrees(bcg.distance_from(xc=g.x,yc=g.y))*60.0*60.0 < tolerance for g in detected_bcg])
            selected = detected_bcg[filter]
            if len(selected) > 0:
                filter = filter[filter]
                bcg_match = bcg_match + 1
                for c in selected:
                    c.parent.bcg_match = True
                    c.parent.match_to = bcg
                    c.parent.separation = degrees(bcg.distance_from(xc=c.x,yc=c.y))*60.0*60.0
                    bcg.matches.append(c.parent)
                    bcg.bcg_match = True
            if len(selected) > 1:
                print 'Double match!'
            
        stringA = 'Found %d direct BCG matches from %d maxBCG clusters'%(bcg_match,len(maxBCG_backup))

        cg_bcg_match = 0
        #STAGE 2: cluster-galaxy maxBCG match - within <tolerance> arcsec

        filter = numpy.array([len(c.matches) == 0 for c in maxBCG])
        maxBCG = maxBCG[filter]

        non_bcg_galaxies = []
        for g in self.cluster_galaxies:
            if g != g.parent.bcg:
                non_bcg_galaxies.append(g)

        for c in self.clusters:
            c.cg_bcg_match = False

        non_bcg_galaxies = numpy.array(non_bcg_galaxies)

        for bcg in maxBCG:
            filter = numpy.array([degrees(bcg.distance_from(xc=g.x,yc=g.y))*60.0*60.0 < tolerance for g in non_bcg_galaxies])
            selected = non_bcg_galaxies[filter]
            if len(selected) > 0:
#                filter = numpy.array([(g.parent.bcg_match == False) or (g.parent.bcg_match and g.parent.match_to != bcg) for g in selected])
                filter = numpy.array([(g.parent.bcg_match == False) for g in selected])
                selected = selected[filter]
            else:
                continue
            print filter    
            if len(selected) > 0:
                filter = filter[filter]
                cg_bcg_match = cg_bcg_match + 1
                for c in selected:
                    if c.parent.bcg_match == False:
                        c.parent.cg_bcg_match = True
                        c.parent.match_to = bcg
                        c.parent.separation = degrees(bcg.distance_from(xc=c.x,yc=c.y))*60.0*60.0
                        bcg.matches.append(c.parent)
                        bcg.cg_bcg_match = True
            if len(selected) > 1:
                print 'Double match!'
            
        stringB = 'Found %d cluster-galaxy:BCG matches from %d maxBCG clusters'%(cg_bcg_match,len(maxBCG_backup))

        #STAGE 2: cluster-galaxy maxBCG match - within <tolerance> arcsec

        maxBCG = maxBCG_backup
        
        #real data
        detector = counter(self.clusters,maxBCG,radii)
        detector = numpy.array(detector,dtype='float')
#        detector = detector/float(len(self.clusters))
        detector = detector/max(detector)

        #random data
        xmin,xmax = self.selection['limits']['xmin'],self.selection['limits']['xmax']
        ymin,ymax = self.selection['limits']['ymin'],self.selection['limits']['ymax']
        dx,dy = xmax-xmin,ymax-ymin
        

        npoints = len(self.clusters)*factor
        all_counts = []
        status = None
        if iterations > 1:
            status = Progress(iterations)

        if r_override == False:
            for i in arange(iterations):
                rx = [xmin+(dx*random.random()) for n in arange(npoints)]
                ry = [ymin+(dy*random.random()) for n in arange(npoints)]
                randoms = []
                for i in xrange(len(rx)):
                    _gal = galaxy(-1,rx[i],ry[i],dummy=True)
                    randoms.append(_gal)
            
                all_counts.append(counter(randoms,maxBCG,radii))
                if status: status.update(0)

            combined_counts = []
            error = []

            normed = [numpy.array(ct,dtype='float')/float(max(ct)) for ct in all_counts]
            for r in radii:
                i = indexOf(radii,r)
                tot = mean([ct[i] for ct in normed])
                e = std([ct[i] for ct in normed])
                combined_counts.append(tot)
                error.append(e)

            
            counts = numpy.array(combined_counts,dtype='float')
            error = numpy.array(error,dtype='float')

            _r = 0.0
            test = [detector[i] > counts[i] for i in xrange(len(radii))]
            loop = iter(test)
            t = loop.next()
            rads = iter(radii)
            while(t == True):
                _r = rads.next()
                t = loop.next()

            _r = _r + dr/2.0

        else: 
            _r = radians(r_override/60.0)

        print 'Using %f-arcmin as the radius cut-off'%(degrees(_r)*60.0)

        filter = numpy.array([len(c.matches) == 0 for c in maxBCG])
        maxBCG = maxBCG[filter]

        for c in self.clusters:
            c.stat_sig = False
            c.null_sig = False
            

        _bcg_matches = numpy.array(self.clusters)[numpy.array([c.bcg_match for c in self.clusters])]
        _cg_bcg_matches = numpy.array(self.clusters)[numpy.array([c.cg_bcg_match for c in self.clusters])]
#        _ss_matches = numpy.array(self.clusters)[numpy.array([c.stat_sig for c in self.clusters])]
#        _null_matches = numpy.array(self.clusters)[numpy.array([c.null_sig for c in self.clusters])]
        _total = list(_bcg_matches)+list(_cg_bcg_matches)
        _total = unique(_total)

        clusters = []
        for c in self.clusters:
            if c not in _total:
                clusters.append(c)

#        filter = numpy.array([(c.bcg_match == False) and (c.cg_bcg_match == False) for c in self.clusters])
#        clusters = numpy.array(self.clusters)[filter]

        clusters = numpy.array(clusters)
        
        stat_sig = 0
        null_counts = 0
        nulls = []

        filter = numpy.array([(True in [c.cg_bcg_match,c.bcg_match,c.stat_sig,c.null_sig]) == False for c in maxBCG])
        maxBCGcl = maxBCG[filter]

        for bcg in maxBCGcl:



            nearby_c = clusters_within_r(bcg,clusters,r=_r)
            within_radius = clusters_within_r(bcg,clusters,r=bcg.extent)
            if len(nearby_c) > 0:
                stat_sig = stat_sig + 1
                if len(nearby_c) > 1:
                    print 'Double ss match!'
                for c in nearby_c:
                    c.stat_sig = True
                    c.match_to = bcg

            elif len(within_radius) > 0:
                _null = []
                for c in within_radius:
                    if (c in nearby_c) == False:
                        _null.append(c)
                        c.null_sig = True
                        c.match_to = bcg
                if len(_null) > 0:
                    null_counts = null_counts + 1
                if len(_null) > 1:
                    print 'Double null match!'
                    
                nulls = nulls + list(_null)
        
        nulls = unique(nulls)
        nulls = len(nulls)

        filter = numpy.array([(c.bcg_match == False) and (c.cg_bcg_match == False) and (c.stat_sig == False) and (c.null_sig == False) for c in self.clusters])
        clusters = numpy.array(self.clusters)[filter]

        _bcg_matches = numpy.array(self.clusters)[numpy.array([c.bcg_match for c in self.clusters])]
        _cg_bcg_matches = numpy.array(self.clusters)[numpy.array([c.cg_bcg_match for c in self.clusters])]
        _ss_matches = numpy.array(self.clusters)[numpy.array([c.stat_sig for c in self.clusters])]
        _null_matches = numpy.array(self.clusters)[numpy.array([c.null_sig for c in self.clusters])]
        _total = list(_bcg_matches)+list(_cg_bcg_matches)+list(_ss_matches)+list(_null_matches)
        _total = unique(_total)

        self.dvorac_matches = _total
        self.bcg_matches = []
        for c in self.dvorac_matches:
            self.bcg_matches.append(c.match_to)
        
            


        clusters = []
        for c in self.clusters:
            if c not in _total:
                clusters.append(c)

        stringC = 'Found %d statistically significant matches to %d maxBCG clusters'%(stat_sig,len(maxBCG_backup))
        stringD = 'Found %d low significance matches to %d maxBCG clusters'%(null_counts,len(maxBCG_backup))
        stringE = 'Found %d clusters that are *not* maxBCG-related at all, from %d detected clusters'%(len(clusters),len(self.clusters))

        print '\n'*2
        print '#### '+stringA
        print '#### '+stringB
        print '#### '+stringC
        print '#### '+stringD
        print '#### '+stringE
        print '\n'*2
        
        results = {}
        results[1] = {'info':stringA,'clusters':_bcg_matches,'maxBCG':numpy.array(unique([c.match_to for c in _bcg_matches]))}
        results[2] = {'info':stringB,'clusters':_cg_bcg_matches,'maxBCG':numpy.array(unique([c.match_to for c in _cg_bcg_matches]))}
        results[3] = {'info':stringC,'clusters':_ss_matches,'maxBCG':numpy.array(unique([c.match_to for c in _ss_matches]))}
        results[4] = {'info':stringD,'clusters':_null_matches,'maxBCG':numpy.array(unique([c.match_to for c in _null_matches]))}





        matches = []
        self.maxBCG = numpy.array(self.maxBCG)
        for c in self.clusters:
            if c.match_to:
                matches.append(c.match_to)

        matches = unique(matches)
        self.matchedBCG_clusters = matches

        unmatched = numpy.array([(c in matches) == False for c in self.maxBCG])
        self.unmatchedBCG_clusters = self.maxBCG[unmatched]

        radii = degrees(radii)*60.0
        error = []
        if plot:
            from detections import plot_clusters
            found = matches
            not_found = self.unmatchedBCG_clusters
            plot_clusters(found,rings=True,centroids=True,colour='!blue')
            plot_clusters(not_found,rings=True,centroids=True,colour='!red',passback=True,limits=self.selection)

            if r_override == False:
                fig = plt.figure()
                if err:
                    plt.errorbar(radii,counts,yerr=error,marker='None',color='r',linestyle='-')
                else:
                    plt.plot(radii,counts,'r-')

                plt.plot(radii,detector,'bo')
                plt.plot([tolerance/60.0,tolerance/60.0],[0,1],'k:')
                
                plt.xlabel('d_theta (arcmins)')
                plt.ylabel('P(theta)')
                plt.ylim(0,1.1)

                

            else:
                counts = []
                error = []

            
        return [radii,results,error,detector]

    def get_ids(self,mode='all',custom=None,verbose=True):
        import numpy
        if mode == 'all' and custom == None:
            gals = numpy.array(self.cluster_galaxies)
        elif mode == 'bcg':
            gals = numpy.array([c.bcg for c in self.clusters])
        elif custom != None:
            gals = numpy.array(custom)


        gid = numpy.array([g.objID for g in gals])
        filter = []
        for id in gid:
            filter.append(id == None)
        filter = numpy.array(filter)

        gals = gals[filter]
        
        if len(gals) == 0:
            return 
        
        gdna = numpy.array([hash('%1.10f %1.10f'%(float(g.x),float(g.y))) for g in gals])

        from dataIO import getData
        dataset = self.selection['io']['dataset']
        catalogue = self.selection['io']['sourcefile']
        obs = self.selection['io'].has_key('observed_data')
        data = getData(catalogue,obs=obs,dataset=dataset,wrap=True)
        x = data['spatial']['x']
        y = data['spatial']['y']
        try:
            id = data['extras']['id']
        except:
            print 'Cannot get objID - aborting'
            return
        target_dna = numpy.array([hash('%1.10f %1.10f'%(float(x[i]),float(y[i]))) for i in xrange(len(x))])
        if verbose:
            print 'Retrieving id matches for objID.....'
        locs = [numpy.where(dna == target_dna)[0] for dna in gdna]

        if verbose:
            print '....done'
        nlocs = numpy.array([len(l) for l in locs])
        filter = (nlocs == 1)
        gals = gals[filter]
        ids = []
        for i in xrange(len(filter)):
            if filter[i]:
                ids.append(id[locs[i][0]])

        for i in xrange(len(gals)):
            gals[i].objID = ids[i]
               

    def determine_arb_boundaries(self,sampling_freq=100):
        from detections import detections as det
        from copy import deepcopy
        from tools import objectify
        from dataIO import Progress
        from clusters import galaxies_within_r
        import numpy
        dummy = det(None,deepcopy(self.selection),dummy=True)
        for colour in dummy.selection['detection']['sequence']:
            dummy.selection['colours'][colour]['width'] = 100.0
        try:
            xx = self.cuts['x']
            xxx = xx[0]
            dummy.cuts = self.cuts
        except:
            dummy.get_input()
#            self.cuts = deepcopy(dummy.cuts)
        dx = dummy.selection['limits']['xmax']-dummy.selection['limits']['xmin']
        dy = dummy.selection['limits']['ymax']-dummy.selection['limits']['ymin']

        dmax = max(dx,dy)
        box_r = dmax/float(sampling_freq)


        x0 = dummy.selection['limits']['xmin']
        patches = []

        dummy.cuts['spatial']['x'] = numpy.array(dummy.cuts['spatial']['x'])
        dummy.cuts['spatial']['y'] = numpy.array(dummy.cuts['spatial']['y'])

        caught = []
        xx = arange(dummy.selection['limits']['xmin'],dummy.selection['limits']['xmax']+box_r/2.0,box_r)
        yy = arange(dummy.selection['limits']['ymin'],dummy.selection['limits']['ymax']+box_r/2.0,box_r)
        
        status = Progress(len(xx)*len(yy))


        for x in xx:
            yy = arange(dummy.selection['limits']['ymin'],dummy.selection['limits']['ymax']+box_r/2.0,box_r)
            this_row = []
            for y in yy:
                xc,yc = x+(box_r/2.0),y+(box_r/2.0)
                _in = galaxies_within_r([xc,yc],dummy.cuts,shape='square',r=box_r)
                if len(_in) >0:
#                    patches.append(plt.Rectangle([xc,yc],box_r,box_r,fc='grey',alpha=0.4,lw=1.0,ec='grey',ls='dashed'))
                    patches.append([[x-(box_r/2.0),y-(box_r/2.0)],box_r,box_r])
                    caught = caught + list(_in)
                    this_row = this_row + list(_in)
                status.update(0)
            if len(this_row) > 0:
                filter = numpy.ones(len(dummy.cuts['spatial']['y']),dtype='bool')
                this_row = numpy.array(this_row)
                filter[this_row] = False
                dummy.cuts['spatial']['x'] = dummy.cuts['spatial']['x'][filter]
                dummy.cuts['spatial']['y'] = dummy.cuts['spatial']['y'][filter]

        self.selection['detection']['boundary'] = {'1':patches}
        return
        

    def apply_darb(self,grid,save=False,status=False,pp=True,ncpu=16,outname=None):
        """\
        det.apply_darb(self,grid,save=False,status=False,pp=True,ncpu=16)
        
        \t Apply a darb grid to the source data this detection belongs to

        grid       :  <string> or <det.grid> containing the boundary data
        save       :  write output file to fits?
        status     :  if pp = False, provide status bar
        pp         :  run checks in parallel
        ncpu       :  #cpu to run un parallel
        """
        dummy = self.clone()
        for colour in dummy.selection['detection']['sequence']:
            dummy.selection['colours'][colour]['width'] = 100.0

        print 'NB: no colour or magnitude restrictons from ORCA - please check limits in .cfg'
        print 'to ensure all input galaxies are appropriate. ORCA will *not* override these!'

        try:
            xx = self.cuts['spatial']['x']
            xxx = xx[0]
#            dummy.cuts = self.cuts
        except:
            try:
                dummy.get_input()
            except:
                print '#detection.apply_darb(): Problem with reading in data'
            
        accept = [False for i in dummy.cuts['spatial']['x']]
        if type(grid) == type(''):
            from dataIO import load
            grid = load(grid)
        dummy.selection['grid'] = grid
        dummy.selection['boundary'] = grid['self']['boundary']
        x = dummy.cuts['spatial']['x']
        y = dummy.cuts['spatial']['y']
        junk = """\
        if status:
            accept = []
            from dataIO import Progress
            status = Progress(len(x))
            for ra,dec in zip(x,y):
                accept.append(dummy.boundary_check(ra,dec))
                status.update(0)
        else:
            accept = [dummy.boundary_check(ra,dec) for ra,dec in zip(x,y)]
        """
        accept = dummy.boundary_check(x,y,status=status,pp=pp,ncpu=ncpu)
        accept = numpy.array(accept)
        filtered = {}
        from tools import objectify,hashify
        obj = numpy.array(objectify(dummy.cuts))
        filtered = obj[accept]
        filtered = hashify(filtered)
        print filtered.keys()
        print dummy.cuts.keys()
        dummy.cuts = filtered
        if save:
            from tools import cuts_to_fits
            if outname == None:
                outname = '%s_barb.fits'%(dummy.selection['io']['dataset'])
            cuts_to_fits(dummy.cuts,outname,det=dummy)
            


        return filtered



    def darb(self,sampling_freq=10,status=True,update_area=True,save=False,area_only=False,healpix=False,nside=128,plot=False,seq_plot=True,full_data=True):
        """\
        detections.darb(self,sampling_freq=10,status=True,update_area=True,save=False,area_only=False) -> [grid,cells,lines,xy,xyid,inout]
        \t DARB : Define ARbitrary Boundary
        Determines the area and boundary of source data. This is important for proper orca analysis of data

       sampling_freq   :  the resolution of the grid - this is the # of subdivisions in the largest of (dRA,dDEC)
       status          :  include status bar?
       update_area     :  update the detection class area, in addition to the area definitions in the <<grid>> data in detection.selection 
       save            :  will save grid data to ./grid_data/<<selection['io']['dataset']>>.grid [will create directory if required - ensure correct permissions!]
       area_only       :  only calculate area, do not determine survey edges
       full_data       :  modify the (dummy) detection class to return *all* galaxies rather than those defined by the .cfg file

        """
        

        from detections import detections as det
        from copy import deepcopy
        from tools import objectify
        from dataIO import Progress
        from clusters import galaxies_within_r,line
        import numpy
        dummy = det(None,deepcopy(self.selection),dummy=True)
        
        def save_grid(detection,grid,empty=False):
            import os
            from dataIO import save
            self = detection
            dataset = self.selection['io']['dataset']
            try:
                stripe = '_stripe'+str(self.selection['detection']['stripe'])
            except:
                stripe = ''
            try:
                if stripe != '':
                    sub_id = 'region'+str(self.selection['detection']['subset_index'])
                else:
                    if self.selection['detection']['subset_index']:
                        sub_id = '_region'+str(self.selection['detection']['subset_index'])
                    else:
                        sub_id = ''
            except:
                sub_id = ''
            if not os.path.isdir('./grid_data'):
                os.mkdir('./grid_data')
                
            if empty:
                grid = {}
                grid['self'] = {'empty':True,'sub_id':self.selection['detection']['subset_index']}
                try:
                    grid['self']['stripe'] = self.selection['detection']['stripe']
                except:
                    pass
            filename = './grid_data/%s%s%s.grid'%(dataset,stripe,sub_id)
            save(grid,filename)

        if full_data:
            print 'Retrieving full data'
            try:
                for colour in dummy.selection['detection']['sequence']:
                    dummy.selection['colours'][colour]['width'] = 100.0
            except:
                pass

        try:
            xx = self.cuts['spatial']['x']
            xxx = xx[0]
            dummy.cuts = self.cuts
        except:
            if 1:
#            try:
                dummy.get_input()
            else:
#            except:
                print '#detection.darb(): Problem with reading in data - cell could be empty?'
                if save:
                    grid = {}
                    save_grid(self,grid,empty=True)
                return None
            



        if len(dummy.cuts['spatial']['x']) == 0:
            print 'Empty data - abandon!!!!'
            if save:
                save_grid(self,grid,empty=True)
            
            return None

        data = deepcopy(dummy.cuts)

        dx = dummy.selection['limits']['xmax']-dummy.selection['limits']['xmin']
        dy = dummy.selection['limits']['ymax']-dummy.selection['limits']['ymin']

        dmax = max(dx,dy)
        box_r = dmax/float(sampling_freq)


        x0 = dummy.selection['limits']['xmin']
        patches = []

        dummy.cuts['spatial']['x'] = numpy.array(dummy.cuts['spatial']['x'])
        dummy.cuts['spatial']['y'] = numpy.array(dummy.cuts['spatial']['y'])

        caught = []
        xx = arange(dummy.selection['limits']['xmin'],dummy.selection['limits']['xmax']+box_r/2.0,box_r)
        yy = arange(dummy.selection['limits']['ymin'],dummy.selection['limits']['ymax']+box_r/2.0,box_r)
        from clusters import vertex,vertex_list

        vertices = {}
        for x in xx:
            for y in yy:
                _v = vertex(0,x,y)
                _v.lines = []
                id = hash('%f %f'%(x,y))
                vertices[id] = {'v':_v,'cells':[],'edge':False}


        if status:
            status = Progress(len(xx)*len(yy))
        lines = {}
        cells = {}
        for x in xx:
            for y in yy:
                hashes = [hash('%f %f'%(x,y)),hash('%f %f'%(x+box_r,y)),hash('%f %f'%(x+box_r,y+box_r)),hash('%f %f'%(x,y+box_r))]
                tester = [h in vertices.keys() for h in hashes]
                verts = []
                if countOf(tester,True) == len(hashes):
                    for h in hashes:
                        try:
                            verts.append(vertices[h]['v'])
                        except:
                            pass

                    cell = vertex_list(verts,None,None,keep_all=True)
                    cell.lines = []
                    id = hash('%f %f'%(x,y))
                    cell.id = id
                    cells[id] = {'cell':cell,'vertices':verts,'full':False}
                    for h in hashes:
                        vertices[h]['cells'].append(cell)


                    l_ids = [hash('%f %f %f %f'%(x,x+box_r,y,y)),hash('%f %f %f %f'%(x+box_r,x+box_r,y,y+box_r)),hash('%f %f %f %f'%(x,x+box_r,y+box_r,y+box_r)),hash('%f %f %f %f'%(x,x,y,y+box_r))]
                    v_l = [[verts[0],verts[1]],[verts[1],verts[2]],[verts[3],verts[2]],[verts[0],verts[3]]]
                    cell_lines = []
                    for lid in l_ids:
                        i = indexOf(l_ids,lid)
                        if lines.has_key(lid):
                            lines[lid].cells.append(cell)
                            cell_lines.append(lines[lid])
#                            if (lines[lid] in lines[lid].vertices[0].lines) == False:
#                                lines[lid].vertices[0].lines.append(lines[lid])
                        else:
                            l = line(v_l[i],cells=[cell])
                            lines[lid]= l
                            cell_lines.append(l)
#                            v_l[i][0].lines.append(l)

                    for l in cell_lines:
                        i = indexOf(cell_lines,l)
                        try:
                            pathway = l.pathway
                            pathway.append(cell_lines[i-3])
                        except:
                            l.pathway = [cell_lines[i-3]]
                                
                    cell.lines = cell_lines

        outside = []
        inside = []
        for x in xx:
            yy = arange(dummy.selection['limits']['ymin'],dummy.selection['limits']['ymax']+box_r/2.0,box_r)
            this_row = []
            for y in yy:
                id = hash('%f %f'%(x,y))
                try:
                    cell = cells[id]
                except:
                    if status:
                        status.update(0)
                    continue
                xc,yc = x+(box_r/2.0),y+(box_r/2.0)
                _in = galaxies_within_r([xc,yc],dummy.cuts,shape='square',r=box_r)
                if len(_in) >0:
#                    patches.append(plt.Rectangle([xc,yc],box_r,box_r,fc='grey',alpha=0.4,lw=1.0,ec='grey',ls='dashed'))
                    patches.append([[x-(box_r/2.0),y-(box_r/2.0)],box_r,box_r])
                    caught = caught + list(_in)
                    this_row = this_row + list(_in)
                    inside.append(cell)
                    cell['cell'].inside = True
                else:
                    outside.append(cell)
                    cell['cell'].inside = False
                if status:
                    status.update(0)
            if len(this_row) > 0:
                filter = numpy.ones(len(dummy.cuts['spatial']['y']),dtype='bool')
                this_row = numpy.array(this_row)
                filter[this_row] = False
                dummy.cuts['spatial']['x'] = dummy.cuts['spatial']['x'][filter]
                dummy.cuts['spatial']['y'] = dummy.cuts['spatial']['y'][filter]

                
        ##### below here is the edge determination. Stop here if you just want area....

        inout = {'in':inside,'out':outside}
        
#        self.darb_area = sum([i['cell'].area() for i in inout['in']])
        self.darb_area = sum([i['cell'].area(sq_element=False) for i in inout['in']])
        estimator = 'DARB'
        if healpix:
            if 1:
#            try:
                from tools import mask
                areas = []
                hp = mask()
                filename = self.selection['io']['sourcefile']
                sub_id = self.selection['grid']['self']['sub_id']
                print filename
                print sub_id
                splitter = '.fits'
                if not splitter in filename and '.fit' in filename:
                    splitter = '.fit'
                    
                

                if '_%s'%(str(sub_id)) not in filename and (sub_id != None):
                    fname = filename.split(splitter)[0]+'_%s.fit'%(sub_id)
#                    fname = filename.split('.')[0]+'_%s.fit'%(sub_id)
                    if os.path.isfile(fname):
                        filename = fname
                    else:
                        print 'not found: %s'%(fname)
                        raise IOError    #this should get us to default to the darb estimate


                if nside == 'auto':
                    nsides = [128,256,512]
                if nside == 'auto_max':
                    nsides = [128,256,512,1024,2048]


                else:
                    nsides = [nside]

                print 'Calculating HEALPIX areas'
                for n in nsides:
                    areas.append(hp.fits_area(filename,x='ra',y='dec',nside=n,units='radians'))
                

#                print self.selection['detection']['detection_area'],update_area,area_only
                
                print 'HEALPIX: %s'%(str(areas))
                print 'DARB   : %s'%(str(self.darb_area))
                if len(areas) == 1:
                    self.darb_area = areas[0]
                    estimator = 'HEALPIX, nside=%d'%(nsides[0])
                else:
                    self.darb_area = numpy.mean(areas)
                    print 'HEALPIX MEAN: %s'%(str(self.darb_area))
                    estimator = 'HEALPIX'
            else:
#            except:
                print 'HEALPix estimate failed - reverting to DARB'

        
        if plot:
            xmax = max([max(c['cell'].vx) for c in cells.values()])+(0.5*box_r)
            xmin = min([min(c['cell'].vx) for c in cells.values()])-(0.5*box_r)
            ymax = max([max(c['cell'].vy) for c in cells.values()])+(0.5*box_r)
            ymin = min([min(c['cell'].vy) for c in cells.values()])-(0.5*box_r)
            xlim = (xmin,xmax)
            ylim = (ymin,ymax)
            col_outside_cell = 'yellow'
            col_inside_cell = 'green'
            col_edge = 'red'
            col_data = 'k'
            spn = 1
            np = 1
            
            if seq_plot:
                from plot import GridPlots
                gap = 0.2
                gp = GridPlots(rows=1, cols=4,horizgap=gap,vertgap=gap,aspect='equal')
                sub_plots = iter(gp.ax)
                fig = gp.fig

#                sps = iter(arange(1,6))
#                spn = sps.next()
#                fig = plt.figure(figsize=plt.figaspect(2/5.0))
                np = 5
          
                
#            if seq_plot:
#                sp = fig.sca(sub_plots.next())
#            else:
#                fig = plt.figure()
#                sp = plt.subplot(1,np,spn,aspect='equal')

#            plt.plot(dummy.cuts['x'],dummy.cuts['y'],col_data+',')
#            for _cell in cells.values():
#                _cell['cell'].plot(colour='k',fill=col_outside_cell,style=':',alpha=0.2)
#            plt.plot(dummy.cuts['x'],dummy.cuts['y'],col_data+',')
#            plt.xlim(xlim)
#            plt.ylim(ylim)
#            sp.axes.get_xaxis().set_ticks([]);sp.axes.get_yaxis().set_ticks([])


            if seq_plot:
                sp = fig.sca(sub_plots.next())
#                spn = sps.next()
            else:
                fig = plt.figure()
                sp = plt.subplot(1,np,spn,aspect='equal')

            for _cell in inside:
                _cell['cell'].plot(colour='k',fill=col_inside_cell,style=':',alpha=0.2)

            for _cell in outside:
                _cell['cell'].plot(colour='k',fill=col_outside_cell,style=':',alpha=0.2)


                plt.plot(dummy.cuts['x'],dummy.cuts['y'],col_data+',')
            plt.xlim(xlim)
            plt.ylim(ylim)
            sp.axes.get_xaxis().set_ticks([]);sp.axes.get_yaxis().set_ticks([])

        if area_only:
            xyid = []
            xy = []

        else:
            for l in lines.values():     
                for v in l.vertices:
                    v.lines.append(l)
            for l in lines.values():     
                if len(l.cells) == 1:
                    cell = l.cells[0]
                    if cell.inside:
                        l.edge = True
                elif len(l.cells) == 2:
                    test = [c.inside for c in l.cells]
                    if countOf(test,True) == 1:
                        l.edge = True
            for l in lines.values():
                l.x = l.vertices[0].x
                l.y = l.vertices[0].y
                if l.edge == True:
                    v1 = l.vertices[1]
                    routes = v1.lines
                    paths = []
                    edge_path1 = []
                    for r in routes:
                        if r != l:
                            paths.append(r)
                            if r.edge == True:
                                edge_path1.append(r)
                    v2 = l.vertices[0]
                    routes2 = v2.lines
                    paths2 = []
                    edge_path2 = []
                    for r in routes2:
                        if r != l:
                            paths2.append(r)
                            if r.edge == True:
                                edge_path2.append(r)
                    l.paths1 = paths
                    l.edge_path1 = edge_path1
                    l.paths2 = paths2
                    l.edge_path2 = edge_path2
            ls = numpy.array(lines.values())
            filter = numpy.array([line.edge for line in ls])
            ls = ls[filter]
            l = ls[1]


            edge_counter = numpy.array([_ll.edge for _ll in lines.values()])
            edges = numpy.array(lines.values())[edge_counter]
            edge_cells = []
            for e in edges:
                edge_cells = edge_cells + e.cells
            n_edges = countOf(edge_counter,True)
            all_regions = []
            all_id = []


            if plot:
                if seq_plot:
                    sp = fig.sca(sub_plots.next())
#                    spn = sps.next()
                else:
                    fig = plt.figure()
                    sp = plt.subplot(1,np,spn,aspect='equal')
                icells = [_c['cell'] for _c in inside]

                for _cell in outside:
                    _cell['cell'].plot(colour='k',fill=col_outside_cell,style=':',alpha=0.2)

                for _cell in inside:
                    _cell['cell'].plot(colour='k',fill=col_inside_cell,style=':',alpha=0.2)
                for _cell in edge_cells:
                    if _cell in icells:
                        _cell.plot(colour=col_edge,fill='None',style=':',alpha=0.5,lw=2.0)

                plt.plot(dummy.cuts['x'],dummy.cuts['y'],col_data+',')
                plt.xlim(xlim)
                plt.ylim(ylim)
                sp.axes.get_xaxis().set_ticks([]);sp.axes.get_yaxis().set_ticks([])



            random_start = True
            do_it = True
            status = Progress(n_edges)
            while(len(all_id) < n_edges):
                print 'Number of edges: %d, Number processed: %d'%(n_edges,len(all_id))
                xyid = []
                xy = []
                filter = numpy.array([lines[key].edge and (key not in all_id) for key in lines.keys()])
                ls = numpy.array(lines.values())[filter]

                j = 1
                if random_start and do_it:
                    jj = arange(1,len(ls)+1)
                    numpy.random.shuffle(jj)
                    j = jj[0]
                    if j == n_edges:
                        j = j -1            ###does this fix a bug?

                    do_it = False
#                    j = 32

                    print 'Using start index %d'%(j)

                    
                l = ls[j]


                xy.append([l.x,l.y])
                xyid.append(l.id)
                status.update(0)


                if plot:
                    counter = 0
                if plot and 0:
                    

                    fig = plt.figure()
                    sp = plt.subplot(1,np,spn,aspect='equal')
                    icells = [_c['cell'] for _c in inside]


                    for _cell in outside:
                        _cell['cell'].plot(colour='k',fill=col_outside_cell,style=':',alpha=0.2)

                    for _cell in inside:
                        _cell['cell'].plot(colour='k',fill=col_inside_cell,style=':',alpha=0.2)
                    for _cell in edge_cells:
                        if _cell in icells:
                            _cell.plot(colour=col_edge,fill='None',style=':',alpha=0.5,lw=2.0)

                    plt.plot(dummy.cuts['x'],dummy.cuts['y'],col_data+',')


#                    plt.plot(dummy.cuts['x'],dummy.cuts['y'],'k,')
                    plt.plot([l.x],[l.y],'b*',ms=20.0)
                    plt.xlim(xlim)
                    plt.ylim(ylim)
                    sp.axes.get_xaxis().set_ticks([]);sp.axes.get_yaxis().set_ticks([])



                while 1:
                    l2 = l.edge_path1[0]
                    if l2.id in xyid:
                        l2 = l.edge_path2[0]
                    if l2.vertices[0] in l.vertices:
                        v = l2.vertices[0]
                    elif l2.vertices[1] in l.vertices:
                        v = l2.vertices[1]
                    _x,_y = v.x,v.y
                    if l2.id in xyid:
                        break

                    xyid.append(l2.id)
                    status.update(0)
                    if not [_x,_y] in xy:
                        xy.append([_x,_y])                    
                        if plot:
                            counter += 1
                            if counter == 8:
                                k = 2

                                if seq_plot:
                                    sp = fig.sca(sub_plots.next())
#                                    spn = sps.next()
                                else:
                                    fig = plt.figure()
                                    sp = plt.subplot(1,np,spn,aspect='equal')
                                icells = [_c['cell'] for _c in inside]

                                for _cell in outside:
                                    _cell['cell'].plot(colour='k',fill=col_outside_cell,style=':',alpha=0.2)


                                for _cell in inside:
                                    _cell['cell'].plot(colour='k',fill=col_inside_cell,style=':',alpha=0.2)
                                for _cell in edge_cells:
                                    if _cell in icells:
                                        _cell.plot(colour=col_edge,fill='None',style=':',alpha=0.5,lw=2.0)



                                plt.plot(dummy.cuts['x'],dummy.cuts['y'],col_data+',')

                                _myxy = numpy.array(xy).transpose()

                                modifier = 1.0
                                _dx = False
                                _dy = False
                                arr_dim = box_r/3.0

                                plt.plot(_myxy[0][:-1],_myxy[1][:-1],'b-',lw=2.0)
                                plt.plot([xy[0][0]],[xy[0][1]],'b*',ms=20.0)

                                m = k+1

                                if _myxy[1][k] != _myxy[1][m]:
                                    _dy = True
                                    if _myxy[1][k] > _myxy[1][m]:
                                        modifier = -1.0

                                if _myxy[0][k] != _myxy[0][m]:
                                    _dx = True
                                    if _myxy[0][k] > _myxy[0][m]:
                                        modifier = -1.0
                                
                                xarr,yarr = _myxy[0][k],_myxy[1][k]

                                if _dy:
                                    yarr += 0.5*box_r
                                    arr = matplotlib.patches.FancyArrow(xarr,yarr,0,modifier*arr_dim,ec='b',fc='b',head_width=arr_dim,head_length=arr_dim,length_includes_head=True)
                                elif _dx:
                                    xarr += 0.5*box_r
                                    arr = matplotlib.patches.FancyArrow(xarr,yarr,modifier*arr_dim,0,ec='b',fc='b',head_width=arr_dim,head_length=arr_dim,length_includes_head=True)
                                    
                                plt.gca().add_patch(arr)


                                _dx = False
                                _dy = False

                                k = -2
                                if _myxy[1][k] != _myxy[1][k-1]:
                                    _dy = True
                                    if _myxy[1][k] < _myxy[1][k-1]:
                                        modifier = -1.0

                                if _myxy[0][k] != _myxy[0][k-1]:
                                    _dx = True
                                    if _myxy[0][k] < _myxy[0][k-1]:
                                        modifier = -1.0

                                
                                xarr,yarr = _myxy[0][k],_myxy[1][k]

                                if _dy:
                                    arr = matplotlib.patches.FancyArrow(xarr,yarr,0,modifier*arr_dim,ec='b',fc='b',head_width=arr_dim,head_length=arr_dim,length_includes_head=True)
                                elif _dx:
                                    arr = matplotlib.patches.FancyArrow(xarr,yarr,modifier*arr_dim,0,ec='b',fc='b',head_width=arr_dim,head_length=arr_dim,length_includes_head=True)
                                    

#                                arr = matplotlib.patches.FancyArrow(xarr,yarr,modifier*arr_dim,0,ec='b',fc='b',head_width=arr_dim,head_length=arr_dim,length_includes_head=True)
                                plt.gca().add_patch(arr)
                               
                                plt.xlim(xlim)
                                plt.ylim(ylim)
                                sp.axes.get_xaxis().set_ticks([]);sp.axes.get_yaxis().set_ticks([])



#                                xarr,yarr = _myxy[0][-1],_myxy[1][-1]
#                                arr = matplotlib.patches.FancyArrow(xarr,yarr,-1.0*arr_dim,0,ec='b',fc='b',head_width=arr_dim,head_length=arr_dim,length_includes_head=True)
#                                plt.gca().add_patch(arr)

#                                xarr,yarr = _myxy[0][2]+(0.5*box_r),_myxy[1][2]+(0.5*box_r)
                                

#                                arr = matplotlib.patches.FancyArrow(_myxy[0][2],_myxy[1][2],0,arr_dim,ec='b',fc='b',head_width=arr_dim,head_length=arr_dim,length_includes_head=True)


#                                arr = plt.Arrow(_myxy[0][2],_myxy[1][2],0,box_r,ec='b',fc='b',width=box_r/2.0,lw=2.0)

#                                plt.show()








                                

                        if xy[0] == xy[-1]:
                            break
                        if len(xy) >= 10*len(ls):
                            break
                    l = l2
                all_regions.append(xy)
                all_id = all_id + xyid



            print 'Finished with a total of %d id'%(len(all_id))
            if plot:
                if seq_plot:
                    sp = fig.sca(sub_plots.next())
#                    spn = sps.next()
                else:
                    fig = plt.figure()
                    sp = plt.subplot(1,np,spn,aspect='equal')

                icells = [_c['cell'] for _c in inside]
                
                for _cell in outside:
                    _cell['cell'].plot(colour='k',fill=col_outside_cell,style=':',alpha=0.2)

                for _cell in inside:
                    _cell['cell'].plot(colour='k',fill=col_inside_cell,style=':',alpha=0.2)
                for _cell in edge_cells:
                    if _cell in icells:
                        _cell.plot(colour=col_edge,fill='None',style=':',alpha=0.5,lw=2.0)
                
                plt.plot(dummy.cuts['x'],dummy.cuts['y'],col_data+',')
                for _xy in all_regions:
#                    _myxy = numpy.array(_xy).transpose()
                    patch = plt.Polygon(_xy,fc='None',alpha=1.0,lw=2.0,ec='blue')
#                    patch = plt.Polygon(_myxy,fc='None',alpha=1.0,lw=2.0,ec='blue')
                    plt.gca().add_patch(patch)
#                    plt.plot(_myxy[0],_myxy[1],'b-',lw=2.0)
#                    conn_x = [_myxy[0][-1],[_myxy[0][0]]]
#                    conn_y = [_myxy[1][-1],[_myxy[1][0]]]
#                    plt.plot(conn_x,conn_y,'b-',lw=2.0)

                    
#                plt.subplots_adjust(wspace=0.0,left=0.02,right=0.98)

                plt.xlim(xlim)
                plt.ylim(ylim)
                sp.axes.get_xaxis().set_ticks([]);sp.axes.get_yaxis().set_ticks([])


                plt.show()

            if len(all_regions) == 1:
                all_regions = all_regions[0]
            self.selection['boundary'] = all_regions
            partitions = []
            if type(self.selection['boundary'][0][0]) == type([]):
                partitions = self.selection['boundary']
                islands = []
            else:
                islands = [self.selection['boundary']]
            holes = []
            from clusters import point_inside_cell
            if len(partitions) > 0:
                print 'Classifying partitions'
                compare = lambda a, b: cmp(len(a),len(b))
                partitions.sort(compare)
                ids = []
                partitions.reverse()

            if len(partitions) > 0:
                    for p in partitions:
                        partitions.reverse()
#                        p_bak = deepcopy(p)
                        for p2 in partitions:
#                            print 'Checking %d against %d'%(len(p),len(p2))
                            if p != p2 and (p2 not in holes):
                                #this can be sped up with a while + break on 1st True....
                                #(done... see .compare for original implementation...)
                                is_hole = False
                                points = iter(p2)
                                while (is_hole == False):
                                    try:
                                        xy = points.next()
                                        is_hole = point_inside_cell(xy[0],xy[1],p,preserve=True)
                                    except StopIteration:
                                        break
                                if is_hole:
                                    holes.append(p2)
                        partitions.reverse()
                    for p in partitions:
                        if not p in holes:
                            #must be an island
                            islands.append(p)
            inout = {'in':inside,'out':outside,'islands':islands,'holes':holes}


        #### cleanup and output:

        from math import pi
        grid = {}
        _this_grid = {}

        if self.selection.has_key('grid') and self.selection['grid'].has_key('self'):
            print 'Recycling existing grid...'
            _this_grid = deepcopy(self.selection['grid']['self'])
        else:
            _this_grid['xmin'] = min(data['spatial']['x'])
            _this_grid['xmax'] = max(data['spatial']['x'])
            _this_grid['ymin'] = min(data['spatial']['y'])
            _this_grid['ymax'] = max(data['spatial']['y'])
            _this_grid['neighbourID'] = []

            try:
                _this_grid['sub_id'] = self.selection['detection']['subset_index']
            except:
                _this_grid['sub_id'] = -1

            try:
                _this_grid['stripe'] = self.selection['detection']['stripe']
            except:
                pass
            n_dna = hash('%1.8f %1.8f'%(xc,yc))
            _this_grid['id'] = n_dna



        junk = """\
        if _this_grid['units'] == 'degrees':
            print 'Converting grid data to radians'
            to_convert = ['xc','yc','xmax','xmin','ymax','ymin']
            for c in to_convert:
                _this_grid[c] = radians(_this_grid[c])

            _this_grid['']
        _this_grid['units'] = 'radians'
        """


        ##### boundary stuff here:

        if area_only == False:

            _this_grid['boundary'] = all_regions
            _this_grid['boundary_class'] = {'island':[],'hole':[]}
            if type(_this_grid['boundary'][0][0]) == type([]):
                for region in _this_grid['boundary']:
                    if region in islands:
                        _this_grid['boundary_class']['island'].append(indexOf(_this_grid['boundary'],region))
                    elif region in holes:
                        _this_grid['boundary_class']['hole'].append(indexOf(_this_grid['boundary'],region))
                    else:
                        print 'Hmm - neither an island nor a hole?'
            self.selection['boundary_class'] = _this_grid['boundary_class']


        ##### spatial/area stuff here
        
        from tools import element_area
        xmin,xmax = _this_grid['xmin'],_this_grid['xmax']
        ymin,ymax = _this_grid['ymin'],_this_grid['ymax']
        xc = xmax-xmin
        yc = ymax-ymin
        dx,dy = xmax-xmin,ymax-ymin
#        _this_grid['area_radians'] = dx*dy
        _this_grid['area_radians'] = element_area(xmin,xmax,ymin,ymax,units='radians')
        


        _this_grid['area_degrees'] = ((360.0*360.0)/(4*pi*pi))*_this_grid['area_radians']       
        _this_grid['no plot'] = True


        grid['self'] = _this_grid
        grid['self']['sampling_frequency'] = sampling_freq
        

        if update_area:
#            print self.selection['detection']['detection_area']

            self.selection['detection']['detection_area_radians'] = self.darb_area
            self.selection['detection']['detection_area'] = self.darb_area*((360.0*360.0)/(4*pi*pi))
            grid['self']['area_radians'] = self.darb_area
            grid['self']['area_degrees'] = self.darb_area*((360.0*360.0)/(4*pi*pi))
            print 'Area for this region (sq deg): %s'%(str(self.selection['detection']['detection_area']))

        self.selection['grid'] = grid
        self.selection['estimator'] = estimator


        if save:
            save_grid(self,grid)

        if not status:
            import gc
#            print 'Arbitrary boundary defined'


            return inout
#        print 'Arbitrary boundary defined'
        return grid,cells,lines,xy,xyid,inout



    def boundary_2d(self,apply=False,plot=False,ngrid=None,verbose=False):
        from tools import dist
        from copy import deepcopy
        stripe_mode = False                                       ## in stripe mode, for apply = True, if a cluster is not allocated a zone, it is added to 'always'                                                                             ## (this accounts for complicated 'always' geometries)
        zones = {}
        zones['always'] = []
        zones['never'] = []
        zones['choose'] = []
        zones['chooser'] = {}
        grid = self.selection['grid']['self']
        self_id = self.selection['grid']['self']['sub_id']
        grids = {}
        for k in self.selection['grid']['self']['neighbourN']:
            id = self.selection['grid'][k]['id']
            if 1:
#            if k != 'self':
                grids[k] = self.selection['grid'][k]
#            else:
#                _id = self.selection['grid']['self']['id']
#                grids[_id] = self.selection['grid']['self']


        if grid['units'] == 'degrees':
            from copy import deepcopy
            _grid = deepcopy(grid)
            keys = ['ddr','dr','r','xc','yc','xmax','xmin','ymax','ymin']
            for key in keys:
                _grid[key] = radians(grid[key])
            grid = _grid



        def dist_check(centre,grids):
            distances = []
            zones = []
            for zone in grids.keys():
                zctr = (grids[zone]['xc'],grids[zone]['yc'])
                distances.append(dist(centre,zctr))
                zones.append(zone)
                
            distances = numpy.array(distances)
            zones = numpy.array([str(z) for z in zones])
#numpy.array(zones)
            
            dlist = deepcopy(distances)
            dlist_srt = deepcopy(dlist)
            dlist_srt.sort()
            thresh = 1e-4
            min_d = min(dlist)
            filter = numpy.array([abs(d-min_d) <= thresh for d in dlist])
            my_zones = zones[filter]
            mz = []
            for z in my_zones:
                try:
                    mz.append(int(z))
                except:
                    mz.append(z)

#            print my_zones
            return list(mz)
            


        def plot_cell(overlap):
            colours = iter(['green','yellow','red','blue'])
#            for g in overlap['grid'].values():
#                rec = plt.Rectangle((g[0],g[1]),g[2],g[3],fc='None',alpha=0.4,lw=2.0,ec='k',ls='dashed')
#                plt.gca().add_patch(rec)
            for g in overlap['never']:
                rec = plt.Rectangle((g[0],g[1]),g[2],g[3],fc='red',alpha=0.4,lw=2.0,ec='red',ls=None)
                plt.gca().add_patch(rec)
            for g in overlap['always']:
                rec = plt.Rectangle((g[0],g[1]),g[2],g[3],fc='green',alpha=0.4,lw=2.0,ec='green',ls=None)
                plt.gca().add_patch(rec)
            for g in overlap['choose']:
                rec = plt.Rectangle((g[0],g[1]),g[2],g[3],fc='yellow',alpha=0.4,lw=2.0,ec='yellow',ls=None)
                plt.gca().add_patch(rec)

            try:
                plt.xlim(overlap['xmin'],overlap['xmin']+overlap['size'])
                plt.ylim(overlap['ymin'],overlap['ymin']+overlap['size'])
            except:
                pass
#            plt.show()
#            raise SystemExit(0)

        dr = grid['dr']
        ddr = 1.0*grid['ddr']
        xc,yc = grid['xc'],grid['yc']
        cell_size = grid['xmax']-grid['xmin']
        
        if plot:
            fig = plt.figure()
            sp = plt.subplot(111,aspect='equal')
            plt.plot([xc],[yc],'kx')
            rec = plt.Rectangle((xc-dr,yc-dr),2.0*dr,2.0*dr,fc='None',alpha=0.0,lw=1.0,ec='k',ls='dotted')
            plt.gca().add_patch(rec)
           

        positions = grid['positions']
        id_lookup = {}
        pos_lookup = {}
        for pos in grid['positions']:
            id_lookup[pos] = grid['neighbourID'][indexOf(grid['positions'],pos)]
            pos_lookup[grid['neighbourID'][indexOf(grid['positions'],pos)]] = pos
        
        positions = []
        ref = arange(-1,1.1,1,dtype='int')
        for r1 in ref:
            for r2 in ref:
                if (r1,r2) != (0,0):
                    positions.append((r1,r2))
                if (r1,r2) in grid['positions']:
                    zones['chooser'][id_lookup[(r1,r2)]] = []
            

        if len(grid['neighbourN']) > 0 and not (grid.has_key('zones') and grid['zones'].has_key('chooser') and grid.has_key('stripe')):              ################### CAN ADD A CAVEAT HERE FOR PRE-COMPUTED STRIPE OVERLAPS? ###############
            xmin = xc-(dr+(1.0/3)*ddr)
            ymin = yc-(dr+(1.0/3)*ddr)
            size = 2*dr+(ddr*2.0/3)
            zones['always'].append((xmin,ymin,size,size))
#            for id in grid['neighbourID']:
#                zones['choose'][id] = []
            dddr = ddr/6.0
            _size = ddr/3.0
            for p in positions:
                sub_grid = {}
                overlap = {}
                overlap['xmin'] = xc-(dr+ddr)
                overlap['ymin'] = yc-(dr+ddr)
                overlap['size'] = 2*dr+(ddr*2.0)
                overlap['never'] = []
                overlap['always'] = []
                overlap['choose'] = []
                overlap['chooser'] = {}
                for k in zones['chooser'].keys():
                    overlap['chooser'][k] = []
                
                if 0 in p:
                    #these are the long sides....
                    if p == (0,1):
                        if p in grid['positions']:
                            overlap['never'].append((xc-dr,yc+(dr+(4*dddr)),dr*2.0,2.0*dddr))
                            overlap['choose'].append((xc-dr,yc+(dr+(2*dddr)),dr*2.0,2.0*dddr))
                            overlap['chooser'][id_lookup[p]].append((xc-dr,yc+(dr+(2*dddr)),dr*2.0,2.0*dddr))
#                            overlap['choose'][p]

                            overlap['always'].append((xc-dr,yc+(dr+(0*dddr)),dr*2.0,2.0*dddr))
                        else:
                            overlap['always'].append((xc-dr,yc+(dr+(0*dddr)),dr*2.0,6.0*dddr))
                    if p == (0,-1):
                        if p in grid['positions']:
                            overlap['never'].append((xc-dr,yc-(dr+(6*dddr)),dr*2.0,2.0*dddr))
                            overlap['choose'].append((xc-dr,yc-(dr+(4*dddr)),dr*2.0,2.0*dddr))
                            overlap['always'].append((xc-dr,yc-(dr+(2*dddr)),dr*2.0,2.0*dddr))
                            overlap['chooser'][id_lookup[p]].append((xc-dr,yc-(dr+(4*dddr)),dr*2.0,2.0*dddr))
                        else:
                            overlap['always'].append((xc-dr,yc-(dr+(6*dddr)),dr*2.0,6.0*dddr))
                    if p == (-1,0):
                        if p in grid['positions']:
                            overlap['never'].append((xc-(dr+(6*dddr)),yc-dr,dddr*2.0,2.0*dr))
                            overlap['choose'].append((xc-(dr+(4*dddr)),yc-dr,dddr*2.0,2.0*dr))
                            overlap['always'].append((xc-(dr+(2*dddr)),yc-dr,dddr*2.0,2.0*dr))
                            overlap['chooser'][id_lookup[p]].append((xc-(dr+(4*dddr)),yc-dr,dddr*2.0,2.0*dr))
                        else:
                            overlap['always'].append((xc-(dr+(6*dddr)),yc-dr,dddr*6.0,2.0*dr))
                    if p == (1,0):
                        if p in grid['positions']:
                            overlap['never'].append((xc+(dr+(4*dddr)),yc-dr,dddr*2.0,2.0*dr))
                            overlap['choose'].append((xc+(dr+(2*dddr)),yc-dr,dddr*2.0,2.0*dr))
                            overlap['always'].append((xc+(dr+(0*dddr)),yc-dr,dddr*2.0,2.0*dr))
                            overlap['chooser'][id_lookup[p]].append((xc+(dr+(2*dddr)),yc-dr,dddr*2.0,2.0*dr))
                        else:
                            overlap['always'].append((xc+(dr+(0*dddr)),yc-dr,dddr*6.0,2.0*dr))

                else:
                    #these are corners
                    x0,y0 = xc+(p[0]*(dr+0.5*ddr)),yc+(p[1]*(dr+0.5*ddr))

                    sub_grid[(0,0)] = (x0-dddr,y0-dddr,_size,_size)
                    no_neighbour = {}
                    centres = {}
                    for r1 in ref:
                        for r2 in ref:
                            _xc = x0+(r1*_size)
                            _yc = y0+(r2*_size)
                            if (r1,r2) == (0,0):
                                centres[(r1,r2)] = (degrees(_xc),degrees(_yc))
                                continue
                            sub_grid[(r1,r2)] = (_xc-dddr,_yc-dddr,_size,_size)
                            centres[(r1,r2)] = (degrees(_xc),degrees(_yc))
                   
                    if 1:
#                    if p in grid['positions']:
                        #changing this to a distance from cell centre criterion,
                        #this should simplify the missing data cell allocation issue...
                        
                        overlap['grid'] = sub_grid
                        overlap['dr'] = dddr
                        n,m = p[0],p[1]
#                        keep = [sub_grid[(-1*n,-1*m)]]
                        reject = []
                        keep = []

                        for key in sub_grid.keys():
                            nearest = dist_check(centres[key],grids)
                            if not 'self' in nearest:
                                reject.append(sub_grid[key])
                            elif len(nearest) == 1 and 'self' in nearest:
                                keep.append(sub_grid[key])
                            elif len(nearest) > 1 and 'self' in nearest:
                                for id in nearest:
                                    if id != 'self':
                                        _id = grids[int(id)]['id']
                                        overlap['choose'].append(sub_grid[key])
                                        overlap['chooser'][_id].append(sub_grid[key])
                                        junk = """\
                                        try:
                                            overlap['chooser'][_id].append(sub_grid[key])
                                        except:
                                            print '\n\nFailed\n\n'
                                            print self.selection['grid']['self']
                                            print '\n\n\n'
                                        """
                            else:
                                print 'wtf? v2'
                                raise NameError('Allocation failed')
                        overlap['never'] = reject
                        overlap['always'] = keep


                        grid_choices = []

                    else:
                        n,m = p[0],p[1]
                        if (n,0) not in grid['positions'] and (0,m) not in grid['positions']:
                            if verbose: print 'missing corner'
                            overlap['always'].append((sub_grid[-1,-1][0],sub_grid[-1,-1][1],ddr,ddr))
                        elif (n,0) not in grid['positions'] and (n,-1*m) not in grid['positions']:
                            overlap['always'].append((sub_grid[(-1,-1*m)][0],sub_grid[(-1,-1*m)][1],ddr,dddr*2))
                            overlap['choose'].append((sub_grid[(-1,0)][0],sub_grid[(-1,0)][1],ddr,dddr*2))
                            overlap['chooser'][id_lookup[(0,m)]].append((sub_grid[(-1,0)][0],sub_grid[(-1,0)][1],ddr,dddr*2))
                            overlap['never'].append((sub_grid[(-1,m)][0],sub_grid[(-1,m)][1],ddr,dddr*2))
                            if verbose: print 'vertical edge, using %s as choice'%(str((0,m)))
                        elif (0,m) not in grid['positions'] and (-1*n,m) not in grid['positions']:
                            overlap['always'].append((sub_grid[(-1*n,-1)][0],sub_grid[(-1*n,-1)][1],dddr*2,ddr))
                            overlap['choose'].append((sub_grid[(0,-1)][0],sub_grid[(0,-1)][1],dddr*2,ddr))
                            overlap['chooser'][id_lookup[(n,0)]].append((sub_grid[(0,-1)][0],sub_grid[(0,-1)][1],dddr*2,ddr))
                            overlap['never'].append((sub_grid[(n,-1)][0],sub_grid[(n,-1)][1],dddr*2,ddr))
                            if verbose: print 'horiz edge, using %s as choice'%(str((n,0)))
                        elif (n,m) not in grid['positions']:
                            overlap['always'].append(sub_grid[p])
                            pass

                        else:
                            if verbose: print 'wtf?'
                            raise SystemExit(0)


                for g in sub_grid.keys():
                    if (sub_grid[g] not in overlap['always']) and (sub_grid[g] not in overlap['never']) and (sub_grid[g] not in overlap['choose']):
                        print 'Should never happen now'
                        raise SystemExit(0)


                if plot:
                    pass
#                    out = plot_cell(overlap)    
                    
                zones['always'] = zones['always'] + list(overlap['always'])
                zones['choose'] = zones['choose'] + list(overlap['choose'])
                if len(list(overlap['choose'])) > 0:
                    for id in overlap['chooser'].keys():
                        if len(overlap['chooser'][id]) > 0:
                            zones['chooser'][id] = zones['chooser'][id] + overlap['chooser'][id]

                zones['never'] = zones['never'] + list(overlap['never'])
                if verbose: print 'done'


            if plot:
                plot_cell(zones)
                if verbose: print (grid['xmin'],grid['xmax']),cell_size,cell_size
                rec = plt.Rectangle((grid['xmin'],grid['ymin']),cell_size,cell_size,fc='None',alpha=0.0,lw=1.0,ec='k',ls='dotted')
                plt.gca().add_patch(rec)
                if ngrid:
                    for id in grid['neighbourID']:

                        nn = ngrid[id]
                        plt.plot([nn['xc']],[nn['yc']],'ro')
                        rec = plt.Rectangle((nn['xmin'],nn['ymin']),cell_size,cell_size,fc='None',alpha=0.0,lw=1.0,ec='r',ls='dotted')
                        plt.gca().add_patch(rec)

                xlim = plt.xlim()
                plt.xlim(xlim[1],xlim[0])
#                plt.show()
            self.zone_layout = zones
#            apply = True

        elif len(grid['neighbourN']) == 0: 
            return

        if (grid.has_key('zones') and grid['zones'].has_key('chooser') and grid.has_key('stripe')):
            stripe_mode = True
            zones = grid['zones']
            if grid['units'] == 'degrees':
                types = grid['zones'].keys()
                for key in types:
                    if key != 'chooser':
                        _tmp_zone = []
                        for zone in grid['zones'][key]:
                            _zone = [radians(p) for p in zone]
                            _tmp_zone.append(_zone)
                        grid['zones'][key] = _tmp_zone
                _chooser = {}
                ids = grid['zones']['chooser'].keys()
                for id in ids:
                    _zone = [radians(p) for p in grid['zones']['chooser'][id]]
                    if type(_zone[0][0]) == numpy.ndarray:
                        print 'Nest fix - change the source grid data to clear this!!!'
                        _zone = _zone[0]
                    _chooser[id] = _zone

                grid['zones']['chooser'] = _chooser
                zones = grid['zones']


        if 1:
            if apply:
                from clusters import vertex,galaxy,vertex_list,point_inside_cell,cluster
                import random
                junk = """\
                clusters = [cluster([],0,dummy=True) for i in arange(100)]
                for cl in clusters:
                    cl.x = grid['xmin']+random.random()*grid['r']
                    cl.y = grid['ymin']+random.random()*grid['r']

                self.clusters = clusters
                """
                never = []
                always = []
                choose = []
                chooser = {}
                sel = {}
                sel['limits'] = {}
                sel['limits']['xmin'],sel['limits']['xmax'] = -10,21
                sel['limits']['ymin'],sel['limits']['ymax'] = 0.7,1.8

                zone_vl = {}
                czones = {}
                czones['always'] = []
                czones['choose'] = []
                czones['never'] = []
                czones['chooser'] = {}
                for cell in zones['never']:
                    x0,y0,dx,dy = cell[0],cell[1],cell[2],cell[3]
                    vertices = [vertex(0,x0,y0),vertex(0,x0+dx,y0),vertex(0,x0+dx,y0+dy),vertex(0,x0,y0+dy)]
                    never.append(vertex_list(vertices,sel,None,keep_all=True))
                    zone_vl['never'] = never
                for cell in zones['always']:
                    x0,y0,dx,dy = cell[0],cell[1],cell[2],cell[3]
                    vertices = [vertex(0,x0,y0),vertex(0,x0+dx,y0),vertex(0,x0+dx,y0+dy),vertex(0,x0,y0+dy)]
                    always.append(vertex_list(vertices,sel,None,keep_all=True))
                    zone_vl['always'] = always
                for id in zones['chooser'].keys():
                    for cell in zones['chooser'][id]:
                        x0,y0,dx,dy = cell[0],cell[1],cell[2],cell[3]
                        vertices = [vertex(0,x0,y0),vertex(0,x0+dx,y0),vertex(0,x0+dx,y0+dy),vertex(0,x0,y0+dy)]
                        _vl = vertex_list(vertices,sel,None,keep_all=True)
                        choose.append(_vl)
                        chooser[_vl] = id
                        
                    zone_vl['choose'] = choose

                tapped = []
                multiple = []
                for key in zone_vl.keys():
                    for zone in zone_vl[key]:
                        for cl in self.clusters:
                            if point_inside_cell(cl.x,cl.y,zone.verts):
                                if cl in tapped:
                                    if cl.myzone in zone_vl[key] == False:
                                        _cl = cl
                                        if verbose: print 'Multiple detection!!!'
                                        plt.plot([cl.x],[cl.y],'kx')
                                        multiple.append(cl)
                                else:
                                    czones[key].append(cl)
                                    tapped.append(cl)
                                    cl.myzone = zone

                if stripe_mode:
                    zoned = []
                    un_zoned = []
                    for cl in self.clusters:
                        try:
                            cl_zone = cl.myzone
                            all_zoned.append(cl)
                        except:
                            un_zoned.append(cl)
                    for cl in un_zoned:
                        czones['always'].append(cl)
                        tapped.append(cl)
                        cl.myzone = zone_vl['always']
                        
                            
                missing = []

                for id in zones['chooser'].keys():
                    czones['chooser'][id] = []
                    

                for cl in self.clusters:
                    if cl.myzone in choose:
                        czones['chooser'][chooser[cl.myzone]].append(cl)

                    if cl not in tapped:
                        if verbose: print 'missing @ %f,%f'%(cl.x,cl.y)
                        plt.plot([cl.x],[cl.y],'k+')

                self.zones = czones

                if plot:
                    colours = {'always':'g','choose':'y','never':'r'}
#                    colours = {'always':'k','choose':'k','never':'k'}
                    for key in czones.keys():
                        if key in colours.keys():
                            for cl in czones[key]:
                                if verbose: print cl.x,cl.y,colours[key]
                                plot_string = '%so'%(colours[key])
                                plt.plot([cl.x],[cl.y],plot_string)
                        

                if verbose: print 'done'



                    



                                                      




    def determine_boundaries(self,thresh=3000.0,wrap=True):
        if self.selection['detection'].has_key('boundary'):
            return
        try:
            data = self.source_data
            xx = data['spatial']['x']
        except:
            from dataIO import getData
            dataset = self.selection['io']['dataset']
            catalogue = self.selection['io']['sourcefile']
            obs = self.selection['io'].has_key('observed_data')
            try:
                wrap = self.selection['io']['wrap']
            except:
                pass
            data = getData(catalogue,obs=obs,dataset=dataset,wrap=wrap)
            self.source_data = data
        x = data['spatial']['x']
        y = data['spatial']['y']
        dx = max(x)-min(x)
        dy = max(y)-min(y)

        constraint = (y > (max(y)-(dy/thresh)))
        _x = x[constraint]
        _y = y[constraint]
        pointA = [min(_x),max(_y)]
        pointB = [max(_x),max(_y)]

        constraint = (y < (min(y)+(dy/thresh)))
        _x = x[constraint]
        _y = y[constraint]
        pointD = [max(_x),max(_y)]
        pointE = [min(_x),max(_y)]

        constraint = (x > (max(x)-(dx/thresh)))
        _x = x[constraint]
        _y = y[constraint]
        pointC = [max(_x),_y[indexOf(_x,max(_x))]]


        constraint = (x < (min(x)+(dx/thresh)))
        _x = x[constraint]
        _y = y[constraint]
        pointF = [min(_x),_y[indexOf(_x,min(_x))]]

        verts = [pointA,pointB,pointC,pointD,pointE,pointF,pointA]
        self.selection['detection']['boundary'] = verts
        return
        
    def to_cluster_zoo(self,target='default',members=False,clobber=False,incomplete_clobber=True,bcg_centred=True):
        """\
        Writes cluster image files to a clusterZoo format ready for analysis
        desintation directory: 
        ./cluster_zoo
        will also write cluster_zoo.fits 
        
        members:
        clobber: overwrite 
        incomplete_clobber: return if *any* files are missing imaging
        bcg_centred: images are centred on the cluster BCG

        """

        import os
        import numpy
        from tools import cat
        from math import degrees
        cwd = os.getcwd()


        if len(self.clusters) == 0:
            return

        if incomplete_clobber:
            test = numpy.array([os.path.isfile(cl.img_name(bcg_mode=bcg_centred)) for cl in self.clusters])
            if False in test:
                dna = numpy.array([cl.dna for cl in self.clusters])
                missing_dna = numpy.where(test == False)[0]
                missing_dna = dna[missing_dna]
                print '\n\nThe following clusterID do not have images:\n'
                print missing_dna
                print '\n\n'
                return

        if target == 'default':
            target = cwd+'/cluster_zoo'
        if os.path.isdir(target) and not clobber:
#            print 'Target %s already exists - bailing'%(target)
            raise SystemExit('Target %s already exists - bailing'%(target))
        elif os.path.isdir(target) and clobber:
            print 'Overwriting target %s '%(target)
            os.system('rm -Rf %s'%(target))
            os.mkdir(target)
        elif not os.path.isdir(target):
            os.mkdir(target)
        fig = plt.figure()
        ax = plt.axes([0,0,1,1], axisbg='white', frameon = True) 
        for cl in self.clusters:
            cl.img_loc = '%s/%d.png'%(target,abs(cl.dna))
            cl.clusterID = cl.dna
            cl.ra = degrees(cl.x)
            cl.dec = degrees(cl.y)
            out = cl.plot(thumb=True,axis=False,grid='AP',passback=True,bcg_centred=bcg_centred);plt.title('')
            plt.subplots_adjust(left=0.0,right=1.00,bottom=0.0,top=1.00)
            plt.savefig(cl.img_loc)
            string = '/usr/bin/convert -fuzz 10% -trim -bordercolor white '+'%s %s'%(cl.img_loc,cl.img_loc)
            print string
            os.system(string)
            plt.clf()
        fits_output = '%s/cluster_zoo.fits'%(target)
        out = cat(self,output=fits_output,explicit=['clusterID','ra','dec','img_loc'])
        
        for cl in self.clusters:
            del cl.img_loc
            del cl.clusterID
            del cl.ra
            del cl.dec

        print 'Done'
        return      
        

    def retrieve_images(self,passback=False,galsize=6,bcgstyle='bs',altscheme=False,offset=0.0,label='',limits='auto',colour=None,fill=None,thumb=False,centroid=False,field_gals=False,naughty=False,factor=0.5,r_max=0.0,bcg_centred=False,grid=True,axis=True,aperture=True,ext=True,aperture_factor=None,rings=False,rgb=None,fieldsize='auto',lw=1.0,grid_clobber=True,send_to='default',mem_conserve=False,copy_to=None,shuffle=False,clobber=False,retry=False,re_render=False):
        import os
        import gc
        thumbs = []
        
        from dataIO import Progress
        status = Progress(len(self.clusters))
        clusters = numpy.array(self.clusters)
        if shuffle:
            numpy.random.shuffle(clusters)
        for c in clusters:
            name = c.img_name()
            if c.has_image(bcg_mode=bcg_centred) and (clobber == False):
                print 'Image already exists for this cluster. To force retrieval, re-sumbit with clobber=True'
                status.update(0)
                continue
            try:
                t = c.retrieve_image(passback=mem_conserve,galsize=galsize,bcgstyle=bcgstyle,altscheme=altscheme,offset=offset,label=label,limits=limits,colour=colour,fill=fill,thumb=thumb,centroid=centroid,field_gals=field_gals,naughty=naughty,factor=factor,r_max=r_max,bcg_centred=bcg_centred,grid=grid,axis=axis,aperture=aperture,ext=ext,aperture_factor=aperture_factor,rings=rings,rgb=rgb,fieldsize=fieldsize,lw=lw,grid_clobber=grid_clobber,send_to=send_to,re_render=re_render)
                if copy_to != None:
                    if copy_to == 'auto':
                        path = '%s.png'%(c.dna)
                        if send_to != 'default':
                            if send_to[-1] != '/':
                                path = '/'+path
                            path = send_to+path
                        plt.title('clusterID=%d'%(c.dna))
                    t.save(full_path=path)
                if mem_conserve == False:
                    thumbs.append(t)
                else:
                    del t
                    ff = gc.collect()
                    
            except:
                print 'Failed for cluster, going for retry'
                if retry == False:
                    continue
                try:
                    t = c.retrieve_image(passback=mem_conserve,galsize=galsize,bcgstyle=bcgstyle,altscheme=altscheme,offset=offset,label=label,limits=limits,colour=colour,fill=fill,thumb=thumb,centroid=centroid,field_gals=field_gals,naughty=naughty,factor=factor,r_max=r_max,bcg_centred=bcg_centred,grid=grid,axis=axis,aperture=aperture,ext=ext,aperture_factor=aperture_factor,rings=rings,rgb=rgb,fieldsize=fieldsize,lw=lw,grid_clobber=grid_clobber,send_to=send_to,re_render=re_render)
                    if copy_to != None:
                        if copy_to == 'auto':
                            path = '%s.png'%(c.dna)
                            if send_to != 'default':
                                if send_to[-1] != '/':
                                    path = '/'+path
                                path = send_to+path
                            plt.title('clusterID=%d'%(c.dna))
                        t.save(full_path=path)
                    if mem_conserve == False:
                        thumbs.append(t)
                    else:
                        del t
                        ff = gc.collect()
                except:
                    pass
            status.update(0)
        return thumbs
                
    def projection_test(self,cluster,n=10,deltaC='auto',keep_cuts=False,job_server=None,alt_mode=False,colour='auto',widen=False,moving_cluster=None):
        from parallel_vertices import deploy_vertices
        from parallel_qhull import qhull
        from parallel_detect import parallel_detect
        from dataIO import sel_edit,pp_run
        import random
        from science import science
        import numpy

        if colour == 'auto':
            try:
                colour = cluster.basis
                magnitude = self.selection['colours'][colour]['magnitude']
            except:
                colour = self.selection['defaults']['colour']
                magnitude = self.selection['colours'][colour]['magnitude']

        if deltaC == 'auto':
            deltaC = cluster.selection['colours'][colour]['width']*3.0
            
        if widen != False:
            self.widen_selection(colour=colour,width=widen,new_input=True)

        sci = science()
        dx = self.selection['limits']['xmax']-self.selection['limits']['xmin']
        dy = self.selection['limits']['ymax']-self.selection['limits']['ymin']

        random_pos = []
        for i in arange(n):
            xr = self.selection['limits']['xmin']+random.random()*dx
            yr = self.selection['limits']['ymin']+random.random()*dy
            random_pos.append([xr,yr])

        if not keep_cuts:
            self.get_input()

        dummy_sel = deepcopy(self.selection)
        det = detections(None,dummy_sel,dummy=True)

        sources = []
        
        cm20_0 = cluster.get_cm20(colour=colour)
        cm20_max = cm20_0 + deltaC
        dcm20 = (cm20_max-cm20_0)/float(n)
        cm20 = arange(cm20_0+dcm20,cm20_max+(dcm20/2.0),dcm20)
        cm20 = [c-cm20_0 for c in cm20]

        det.get_input()
        
        target = cluster
        duplicate = True
        insert = False
        shunt = True
        dummy_array = numpy.array([])

        if moving_cluster != None:

            if type(moving_cluster) == type([]) or type(moving_cluster) == type(dummy_array):
                target_cluster_list = []
                for cl in moving_cluster:
                    cl.parent = None
                    target = deepcopy(cl)
                    target.move_to(x=cluster.x,y=cluster.y)
                    rmax = cluster.max_radius()
                    rx = random.random()-0.5
                    ry = random.random()-0.5
                    target.move_to(x=cluster.x+(rx*rmax),y=cluster.y+(ry*rmax))
                    target_cluster_list.append(target)
                    

            else:
                target = deepcopy(moving_cluster)
                target.parent = None
                target.move_to(x=cluster.x,y=cluster.y)
                rmax = cluster.max_radius()
                rx = random.random()-0.5
                ry = random.random()-0.5
                target.move_to(x=cluster.x+(rx*rmax),y=cluster.y+(ry*rmax))
                target_cluster_list = [target]


            duplicate = False
            insert = True
            shunt = False
                
            
        else:
            moving_cluster = target
            target = deepcopy(cluster)
            target.parent = None
            target.move_to(x=cluster.x+radians(1.0/60.0),y=cluster.y+radians(1.0/60.0))
#            cl_dict = hashify(target.galaxies)
#            merge_data(det.cuts,cl_dict)
            duplicate = False
            insert = True
            shunt = False
            target_cluster_list = [target]
            
        target_clusters = []



#        target.parent = None
        print 'Generating cluster list'
        source = deepcopy(det)
        source.set_tracer(cluster,traceid='origin')
        source.cuts['extras']['tracer'] = numpy.array([False for f in source.cuts['spatial']['x']])

        tcl = []
        for tc in target_cluster_list:
            tc.parent = None
            target_clusters = []
            for cm in cm20:
                _tc = deepcopy(tc)
                for g in _tc.galaxies:
                    g.colours[colour] = g.colours[colour] + cm
                    g.extras['tracer'] = True
                    g.extras['origin'] = False
                target_clusters.append(_tc)
            tcl.append(target_clusters)
        target_cluster_list = tcl
        if 1:
            pass

        else:
            if 1:
                source.cuts = self.project_to(target,0.0,cm,do_detection=False,return_source=True,keep_cuts=True,tracer=True,colour=colour,displacement='colour',shunt=shunt,duplicate=duplicate,insert=insert)
                source.set_tracer(cluster,traceid='origin')            
                source.get_input(data=source.cuts)
            sources.append(source)
            print 'Done 1'
        dependent_functions = (deploy_vertices,qhull,parallel_detect,sel_edit,sci.richness)
        modules = ('science','parallel_detect','clusters','numpy','parallel_qhull','parallel_vertices','vperc','time','detections')
        run_function = 'parallel_detect'
        args = [{'source[i].cuts':None},{'source[i].selection':None}]            
        ncpus = 16
        print 'Going for detections....'
        completed,job_server = pp_run([source],run_function,args,dependent_functions,modules,ncpus=ncpus,idle_delay=360.0,verbose=False,server=job_server)

        fracs = []

        for target_clusters in target_cluster_list:
            sources = []
            dummy_sel = deepcopy(self.selection)
            det = detections(None,dummy_sel,dummy=True)

            for i in xrange(len(cm20)):
                det = detections(None,dummy_sel,dummy=True)
                clusters = completed[0]+[target_clusters[i]]
                det.addclusters(clusters,verbose=False)
                sources.append(det)
            
            alpha = []
            beta = []

            separated = False
            dets = iter(sources)
            cmseps = iter(cm20)
            from dataIO import Progress
            print 'Checking cluster with %d members vs mobile cluster with %d members'%(cluster.ngalaxies,target_clusters[0].ngalaxies)
            status = Progress(len(sources))

            while(separated == False):
                try:
                    det = dets.next()
                    dc = cmseps.next()
                except:
                    print 'Did not resolve out in time - increase deltaC'
                    dc = self.selection['colours'][colour]['width']*2.0
                    break
                is_sep = False
            
                det.clusters = numpy.array(det.clusters)
                origin = numpy.array([True in [g.extras['origin'] for g in c.galaxies] for c in det.clusters])
                tracer = numpy.array([True in [g.extras['tracer'] for g in c.galaxies] for c in det.clusters])
                origin = det.clusters[origin]
                tracer = det.clusters[tracer]
                if len(origin) > 0:
                    origin = origin[0]
                    tracer = tracer[0]
                    match = origin.matches([tracer],cm20_thresh=self.selection['colours'][colour]['width']*2.0)
                    if len(match) == 0:
                        separated = True
                status.update(0)
            frac = dc/(self.selection['colours'][colour]['width']*2.0)
            print 'Resolved out @ dc=%f,which is %f percent of the filterwidth'%(dc,frac*100.0)
            fracs.append(frac)

        return [fracs,dc,sources],job_server


    def projection_test2(self,cluster,n=10,deltaC='auto',keep_cuts=False,job_server=None,alt_mode=False,colour='auto',widen=False,moving_cluster=None):
        from parallel_vertices import deploy_vertices
        from parallel_qhull import qhull
        from parallel_detect import parallel_detect
        from dataIO import sel_edit,pp_run
        import random
        from science import science
        import numpy

        if colour == 'auto':
            try:
                colour = cluster.basis
                magnitude = self.selection['colours'][colour]['magnitude']
            except:
                colour = self.selection['defaults']['colour']
                magnitude = self.selection['colours'][colour]['magnitude']

        if deltaC == 'auto':
            deltaC = cluster.selection['colours'][colour]['width']*3.0
            
        if widen != False:
            self.widen_selection(colour=colour,width=widen,new_input=True)

        sci = science()
        dx = self.selection['limits']['xmax']-self.selection['limits']['xmin']
        dy = self.selection['limits']['ymax']-self.selection['limits']['ymin']

        random_pos = []
        for i in arange(n):
            xr = self.selection['limits']['xmin']+random.random()*dx
            yr = self.selection['limits']['ymin']+random.random()*dy
            random_pos.append([xr,yr])

        if not keep_cuts:
            self.get_input()

        dummy_sel = deepcopy(self.selection)
        det = detections(None,dummy_sel,dummy=True)

        sources = []
        
        cm20_0 = cluster.get_cm20(colour=colour)
        cm20_max = cm20_0 + deltaC
        dcm20 = (cm20_max-cm20_0)/float(n)
        cm20 = arange(cm20_0+dcm20,cm20_max+(dcm20/2.0),dcm20)
        cm20 = [c-cm20_0 for c in cm20]

        
        target = cluster
        duplicate = True
        insert = False
        shunt = True
        if moving_cluster:
            target = moving_cluster
            target.move_to(x=cluster.x,y=cluster.y)

            if colour != self.selection['detection']['sequence'][-1]:
                #need to change the cluster's next colour sequence so it remains detectable...
                _next_col = self.selection['detection']['sequence'][indexOf(self.selection['detection']['sequence'],colour)+1]
                _next_mag = self.selection['colours'][_next_col]['magnitude']
                _mags = iter(self.cuts['magnitudes'][_next_mag])
                _cols = iter(self.cuts['colours'][_next_col])
                for g in target.galaxies:
                    g.magnitudes[_next_mag] = _mags.next()
                    g.colours[_next_col] = _cols.next()

            duplicate = False
            insert = True
            shunt = False
        for cm in cm20:
            source = deepcopy(det)
            source.cuts = self.project_to(target,0.0,cm,do_detection=False,return_source=True,keep_cuts=True,tracer=True,colour=colour,displacement='colour',shunt=shunt,duplicate=duplicate,insert=insert)
            source.set_tracer(cluster,traceid='origin')
            source.get_input(data=source.cuts)
            sources.append(source)
    
        dependent_functions = (deploy_vertices,qhull,parallel_detect,sel_edit,sci.richness)
        modules = ('science','parallel_detect','clusters','numpy','parallel_qhull','parallel_vertices','vperc','time','detections')
        run_function = 'parallel_detect'
        args = [{'source[i].cuts':None},{'source[i].selection':None}]            
        ncpus = 16
        print 'Going for detections....'
        completed,job_server = pp_run(sources,run_function,args,dependent_functions,modules,ncpus=ncpus,idle_delay=360.0,verbose=False,server=job_server)

        for i in xrange(len(sources)):
            sources[i].addclusters(completed[i])
            
        alpha = []
        beta = []
        
        separated = False
        dets = iter(sources)
        cmseps = iter(cm20)
        while(separated == False):
            try:
                det = dets.next()
                dc = cmseps.next()
            except:
                print 'Did not resolve out in time - increase deltaC'
                dc = -99.0
                break
            is_sep = False
            for c in det.clusters:
                if (True in [g.extras['origin'] for g in c.galaxies]):
                    if (True in [g.extras['tracer'] for g in c.galaxies]) == False:
                        is_sep = True

            separated = is_sep

            
        frac = dc/(self.selection['colours'][colour]['width']*2.0)
        print 'Resolved out @ dc=%f,which is %f percent of the filterwidth'%(dc,frac*100.0)
        return [frac,dc,sources]












    def robustness(self,trainer,per_bin=10,plot=False,dfrac=0.1,mode='median'):
        from detections import detections
        from clusters import cluster
        import numpy
        from dataIO import sel_edit2

        tdna = numpy.array([hash('%1.8f %1.8f'%(float(g.x),float(g.y))) for g in trainer.galaxies])

        self.selection['limits']['magnitudes']['r']['limit'] = 23.5

        try:
            data = self.source_data
            xx = data['spatial']['x']
        except:
            from dataIO import getData
            dataset = self.selection['io']['dataset']
            catalogue = self.selection['io']['sourcefile']
            obs = self.selection['io'].has_key('observed_data')
#            data = getData(catalogue,obs=obs,dataset=dataset,wrap=False)
#            self.source_data = data

        try:
            del self.selection['detection']['subset_index']
        except:
            pass

        trainer.parent.selection = self.selection


        #make any selection changes here if required...

        tr_index,field,shared_dna,t_gals = trainer.sel_completeness(precision=8,verbose=True)

        data = trainer.parent.cuts

#        source,js = sel_edit2(self.selection,'gr_width',0.3,perfect_clusters=False,pp=False)
        
        x = data['spatial']['x']
        y = data['spatial']['y']
        source_dna = numpy.array([hash('%1.8f %1.8f'%(float(x[i]),float(y[i]))) for i in xrange(len(x))])
        junk = """\
        field = []
        tr_index = []
        _in_source = []
        for i in xrange(len(source_dna)):
            if source_dna[i] in tdna:
                tr_index.append(i)
                _in_source.append(source_dna[i])
            else:
                field.append(i)

        field = numpy.array(field)
        tr_index = numpy.array(tr_index)

        t_gals = []
        for d in in_source:
            _ii = indexOf(tdna,d)
            t_gals.append(trainer.galaxies[_ii])

        """
        field_cuts = {}
        for key in data.keys():
            if type(data[key]) == type({}):
                field_cuts[key] = {}
                for subkey in data[key].keys():
                    try:
                        field_cuts[key][subkey] = data[key][subkey][field]
                    except:
                        field_cuts[key][subkey] = None
            else:
                field_cuts[key] = data[key][field]
        


        

        dummy_sel = deepcopy(self.selection)
        det = detections(None,dummy_sel,dummy=True)

        test_cluster = cluster(t_gals,0,real=True)
        test_det = detections(None,dummy_sel,dummy=True)
        test_det.addclusters([test_cluster])
        gals = []
        frac_rej = []


        string = """\

rescale how figure.py is doing by increasing width of g-r. Can we actually detect
the all the galaxies in that selection? Increase until we do. Is the .cmr() not doing
the right thing? r_mag max ~23.3, so change the limiting magnitude too. Then re-run the 
4x4 plot. Can we re-write the above as a function: for a given selection, what galaxies
in a cluster can you even see (need to worry about hashing precision...)


        """
        bcg_dna = hash('%1.8f %1.8f'%(trainer.bcg.x,trainer.bcg.y))
        bcg_index = indexOf(source_dna,bcg_dna)
        

        from parallel_vertices import deploy_vertices
        from parallel_qhull import qhull
        from parallel_detect import parallel_detect
        from dataIO import sel_edit,pp_run
        from science import science
        sci = science()
        dependent_functions = (deploy_vertices,qhull,parallel_detect,sel_edit,sci.richness)
        modules = ('science','parallel_detect','clusters','numpy','parallel_qhull','parallel_vertices','vperc','time','detections')
        run_function = 'parallel_detect'
        args = [{'source[i].cuts':None},{'source[i].selection':None}]            
        job_server = None
        ncpus = 16

        completeness = []
        purity = []

        ejection_fraction = arange(0.0,0.96,dfrac)
#        ejection_fraction = [0.4]
        for f in ejection_fraction:
            sources = []
            realisations = arange(0,per_bin)
            f_dets,f_fields = [],[]
            ngals = int((1.0-f)*float(len(tr_index)))
            frac_rej.append(f)
            training_detections = []
            for case in realisations:
                src = deepcopy(det)
                _this_cluster = deepcopy(tr_index)
                numpy.random.shuffle(_this_cluster)
                sequence = iter(_this_cluster)
                indices = []
                indices.append(bcg_index)
                while(len(indices) != ngals):
                    try:
                        j = sequence.next()
                    except:
                        print 'Have run out of available indices'
                        print 'Using fraction %f. Need %d gals from trainer index length %d'%(f,ngals,len(_this_cluster))
                        print 'Currently have %d in list'%(len(indices))
                        j = sequence.next()
                    if j not in indices:
                        indices.append(j)
                indices = numpy.array(indices)
    
                linkage = [indexOf(tr_index,_t) for _t in indices]

                selected_trainer_galaxies = t_gals[linkage]
                this_trainer = cluster(selected_trainer_galaxies,0,real=True)
                _det = deepcopy(test_det)
                _det.addclusters([this_trainer])
                training_detections.append(_det)
                cuts = {}
                for key in data.keys():
                    if type(data[key]) == type({}):
                        cuts[key] = {}
                        for subkey in data[key].keys():
                            try:
                                cuts[key][subkey] = numpy.array(list(data[key][subkey][indices])+list(field_cuts[key][subkey]))
                            except:
                                cuts[key][subkey] = None
                    else:
                        cuts[key] = numpy.array(list(data[key][indices])+list(field_cuts[key]))
                src.cuts = cuts        
                sources.append(src)
            completed,job_server = pp_run(sources,run_function,args,dependent_functions,modules,ncpus=ncpus,idle_delay=360.0,verbose=False,server=job_server)
            for case in completed:
                _j = indexOf(completed,case)
                src.addclusters(case)
                f_detected,f_field,clust = src.performance(training_detections[_j],plot=False)
                f_dets.append(f_detected)
                f_fields.append(float(len(test_cluster.galaxies))/float(len(clust.galaxies)))
#                f_fields.append(1.0-f_field)
            completeness.append(f_dets)
            purity.append(f_fields)
        
        c,cerr,p,perr = [],[],[],[]


        if mode == 'mean':
            from numpy import mean as median

        from numpy import median,std
        for i in xrange(len(completeness)):
            c.append(median(completeness[i]))
            cerr.append(std(completeness[i]))
            p.append(median(purity[i]))
            perr.append(std(purity[i]))

        output = {}
        output['rejection fraction'] = frac_rej
        output['c'] = c
        output['cerr'] = cerr
        output['p'] = p
        output['perr'] = perr

        if plot:
            plt.errorbar(output['rejection fraction'],output['c'],yerr=output['cerr'],marker='o',color='b',linestyle='-',label='Completeness')
            plt.errorbar(output['rejection fraction'],output['p'],yerr=output['perr'],marker='o',color='r',linestyle='-',label='Purity')
            plt.legend(loc=1)
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.xlabel('Rejection fraction')

        return output
        


    def implant_cluster(self,cluster,x='auto',y='auto',do_detection=True,return_source=False,keep_cuts=False,mode='cuts'):
        from parallel_detect import parallel_detect
        from tools import objectify,hashify

        input = hashify(cluster.galaxies)
        if not keep_cuts:
            self.get_input()

        if x != 'auto':
            _x = cluster.x-x
            input['spatial']['x'] = input['spatial']['x'] - _x

        if y != 'auto':
            _y = cluster.y-y
            input['spatial']['y'] = input['spatial']['y'] - _y


        if mode == 'cuts':
            dummy_sel = deepcopy(self.selection)
            det = detections(None,dummy_sel,dummy=True)
            det.cuts = deepcopy(self.cuts)
            source = det.cuts

        elif mode == 'data':
            source = deepcopy(self.source_data)

        for key in source.keys():
            if input.has_key(key):
                if type(source[key]) == type({}):
                    for subkey in source[key].keys():
                        if input[key].has_key(subkey):
                            try:
                                source[key][subkey] = numpy.array(list(source[key][subkey])+list(input[key][subkey]))
                            except:
                                source[key][subkey] = None
                                
                else:
                    source[key] = numpy.array(list(source[key])+list(input[key]))

        if do_detection and mode == 'cuts':
            output = parallel_detect(det.cuts,det.selection)
            det.addclusters(output)
            return det
        if return_source:
            return source
            

    def displacement_test(self,cluster,n=10,keep_cuts=False,job_server=None,alt_mode=False,buffer_factor=1.5,edge_mode=False):
        from parallel_vertices import deploy_vertices
        from parallel_qhull import qhull
        from parallel_detect import parallel_detect
        from dataIO import sel_edit,pp_run
        import random
        from science import science
        sci = science()
        dx = self.selection['limits']['xmax']-self.selection['limits']['xmin']
        dy = self.selection['limits']['ymax']-self.selection['limits']['ymin']

        r_max = cluster.max_radius()
        x_range = (self.selection['limits']['xmin']+(buffer_factor*r_max),self.selection['limits']['xmax']-(buffer_factor*r_max))
        y_range = (self.selection['limits']['ymin']+(buffer_factor*r_max),self.selection['limits']['ymax']-(buffer_factor*r_max))

        
        
        random_pos = []
        if edge_mode:
            while len(random_pos) != n:
                #choose region of survey to move to
#                corners = [1,3,5,7]
                corners = []
                modes = [1,2,3,4]
                random.shuffle(modes)
                mode = modes[0]
                range_x = [self.selection['limits']['xmin']+r_max,self.selection['limits']['xmax']-r_max]
                range_y = [self.selection['limits']['ymin']+r_max,self.selection['limits']['ymax']-r_max]
                cdx = range_x[1]-range_x[0]
                cdy = range_y[1]-range_y[0]
                if 1:
                    if mode == 1:
                        xr = range_x[0]+(random.random()*cdx)
                        yr = range_y[1]-((buffer_factor-1.0)*random.random())*r_max
                    elif mode == 2:
                        yr = range_y[0]+(random.random()*cdy)
                        xr = range_x[1]-((buffer_factor-1.0)*random.random())*r_max
                    elif mode == 3:
                        xr = range_x[0]+(random.random()*cdx)
                        yr = range_y[0]+((buffer_factor-1.0)*random.random())*r_max
                    elif mode == 4:
                        yr = range_y[0]+(random.random()*cdy)
                        xr = range_x[0]+((buffer_factor-1.0)*random.random())*r_max
                random_pos.append([xr,yr])
            junk = """\
            xy = (self.selection['limits']['xmin'],self.selection['limits']['ymin'])
            rec = plt.Rectangle(xy,dx,dy,fc='green',alpha=0.4,lw=2.0,ec='green',ls='dashed')
            plt.gca().add_patch(rec)
            for p in random_pos:
                plt.plot([p[0]],[p[1]],'ro')
            plt.show()
            raise SystemExit(0)
            """        
                    


        while len(random_pos) != n:
            xr = self.selection['limits']['xmin']+random.random()*dx
            yr = self.selection['limits']['ymin']+random.random()*dy
            while(xr < x_range[0]) or (xr > x_range[1]):
                xr = self.selection['limits']['xmin']+random.random()*dx
            while(yr < y_range[0]) or (yr > y_range[1]):
                yr = self.selection['limits']['ymin']+random.random()*dy

            random_pos.append([xr,yr])

        print 'Generated random positions'
            

        if not keep_cuts:
            self.get_input()

        dummy_sel = deepcopy(self.selection)
        det = detections(None,dummy_sel,dummy=True)

        sources = []
        for r in random_pos:
            source = deepcopy(det)
            source.cuts = self.project_to(cluster,r[0],r[1],do_detection=False,return_source=True,keep_cuts=True,tracer=True)
            sources.append(source)
        
        dependent_functions = (deploy_vertices,qhull,parallel_detect,sel_edit,sci.richness)
        modules = ('science','parallel_detect','clusters','numpy','parallel_qhull','parallel_vertices','vperc','time','detections')
        run_function = 'parallel_detect'
        args = [{'source[i].cuts':None},{'source[i].selection':None}]            
        ncpus = 16
        print 'Going for detections....'
        completed,job_server = pp_run(sources,run_function,args,dependent_functions,modules,ncpus=ncpus,idle_delay=360.0,verbose=False,server=job_server)

        for i in xrange(len(sources)):
            sources[i].addclusters(completed[i])

        completeness = []
        gal_frac = []
        purity = []
        from clusters import clusters_within_r,galaxies_within_r,galaxy
        nparents = []
        member_fraction = []
        for src in sources:
            j = indexOf(sources,src)
            trace = numpy.array([g.extras['tracer'] for g in src.cluster_galaxies])
            galaxies = numpy.array(src.cluster_galaxies)
            traced = galaxies[trace]
#            if len(traced) == 1.1:

            parents = unique([g.parent for g in traced])
            if len(parents) > 0:
                frac = [countOf([g.extras['tracer'] for g in p.galaxies],True)/float(p.ngalaxies) for p in parents]
                j = indexOf(frac,max(frac))
                _p = parents[j]
                mf = float(float(len(p.galaxies))/countOf([g.extras['tracer'] for g in _p.galaxies],True))
                
#                member_fraction.append(mean(frac))
                member_fraction.append(mf)

            nparents.append(len(parents))




            if 1:
                cgals = galaxies_within_r(galaxy(0,random_pos[j][0],random_pos[j][1],dummy=True),src.clusters,r=radians(10.0/60.0))
#                parents = unique([c.parent for c in cgals])

                completeness.append(float(len(traced))/float(len(cluster.galaxies)))
                if len(parents) == 1:                
                    _p = parents[0]
                purity.append(float(len(cluster.galaxies))/float(len(_p.galaxies)))



#                if len(traced) >= int(0.5*len(cluster.galaxies)):
#                if len(clusters_within_r(galaxy(0,random_pos[j][0],random_pos[j][1],dummy=True),src.clusters,r=radians(30.0/60.0))) >0:
#                    completeness.append(1.0)
#                elif len(traced) > 1:
#                if len(clusters_within_r(galaxy(0,random_pos[j][0],random_pos[j][1],dummy=True),src.clusters,r=radians(30.0/60.0))) >0:
#                    completeness.append(0.1)
#                else:
#                    completeness.append(0.0)
#                purity.append(0.0)
                gal_frac.append(0.0)
                continue

          

            if len(parents) == 1.141:
                parent = parents[0]
                if alt_mode:
                    if len(clusters_within_r(galaxy(0,random_pos[j][0],random_pos[j][1],dummy=True),r=radians(3.0/60.0))) >0:
#                    if len(traced) >= 1:
#                    if len(traced) >= 5 or (len(trace) == 5 and len(traced) >= 3):
                        completeness.append(1.0)
                else:
                    completeness.append(float(len(traced))/float(len(cluster.galaxies)))

                traced = numpy.array(parent.galaxies)
                filter = numpy.array([g.extras['tracer'] for g in parent.galaxies])
                traced = traced[filter]

                n_fd = len(parent.galaxies)-len(traced)
                n_det = len(parent.galaxies)
                purity.append(1.0-(float(n_fd)/float(n_det)))
                gal_frac.append(float(len(traced))/float(n_det))

            elif len(parents) == 2.141:
                counts = []
#                for p in parents:
#                    trace_members = [g.extras['tracer'] for g in p.galaxies]
#                    counts.append(countOf(trace_members,True))
#                parent = parents[indexOf(counts,max(counts))]

                if alt_mode:
                    if len(clusters_within_r(galaxy(0,random_pos[j][0],random_pos[j][1],dummy=True),r=radians(3.0/60.0))) >0:
#                    if len(traced) >= 1:
#                    if len(traced) >= 5 or (len(trace) == 5 and len(traced) >= 3):
                        completeness.append(1.0)
                else:
                    completeness.append(float(len(traced))/float(len(cluster.galaxies)))
   #                completeness.append(float(max(counts))/float(len(cluster.galaxies)))

                traced = numpy.array(parent.galaxies)
                filter = numpy.array([g.extras['tracer'] for g in parent.galaxies])
                traced = traced[filter]

                n_fd = len(parent.galaxies)-len(traced)
                n_det = len(parent.galaxies)
                purity.append(1.0-(float(n_fd)/float(n_det)))
                gal_frac.append(float(len(traced))/float(n_det))
        if len(member_fraction) == 0:
            member_fraction = [0.0]
            
        return [completeness,purity,gal_frac,job_server,member_fraction,nparents]


    def set_tracer(self,cluster,traceid='tracer'):
        try:
            a = len(self.cuts['spatial']['x'])
        except:
            self.get_input()
            
        source = self.cuts
        datax = source['spatial']['x']
        datay = source['spatial']['y']

        cg_dna = [hash('%1.8f %1.8f'%(g.x,g.y)) for g in cluster.galaxies]
        source_dna = [hash('%1.8f %1.8f'%(datax[i],datay[i])) for i in xrange(len(datax))]

        to_tag = []
        in_source = []
        to_move = []
        for i in xrange(len(cg_dna)):
            if cg_dna[i] in source_dna:
                to_move.append(indexOf(source_dna,cg_dna[i]))
                in_source.append(cluster.galaxies[i])

        if len(to_move) != len(cluster.galaxies):
            print 'Warning - incomplete detection of cluster galaxies in source catalogue!'
        
        if not source.has_key('extras'):
            source['extras'] = {}
        data_range = arange(len(source['spatial']['x']))
        trace = numpy.array([_d in to_move for _d in data_range])
        source['extras'][traceid] = trace
            





    def project_to(self,cluster,x,y,do_detection=True,return_source=False,keep_cuts=False,mode='cuts',tracer=False,displacement='spatial',colour='auto',shunt=False,duplicate=False,insert=False):
        from parallel_vertices import deploy_vertices
        from parallel_qhull import qhull
        from parallel_detect import parallel_detect
        from dataIO import sel_edit,pp_run
        from science import science
        sci = science()
        dependent_functions = (deploy_vertices,qhull,parallel_detect,sel_edit,sci.richness)
        modules = ('science','parallel_detect','clusters','numpy','parallel_qhull','parallel_vertices','vperc','time','detections')
        run_function = 'parallel_detect'
        args = [{'source[i].cuts':None},{'source[i].selection':None}]            
        job_server = None
        ncpus = 16

        from tools import objectify,hashify
        from clusters import galaxy
        from clusters import galaxies_within_r
        nclusters = len(self.clusters)
        if not keep_cuts:
            self.get_input()

            
        if mode == 'cuts':
            dummy_sel = deepcopy(self.selection)
            det = detections(None,dummy_sel,dummy=True)
            det.cuts = deepcopy(self.cuts)
            source = det.cuts

        elif mode == 'data':
            source = deepcopy(self.source_data)

        
        if insert:
            from tools import merge_data
            cl_dict = hashify(cluster.galaxies)
            merge_data(source,cl_dict)


        dummy_sel = deepcopy(self.selection)
        det = detections(None,dummy_sel,dummy=True)
        det.cuts = deepcopy(self.cuts)

        datax = source['spatial']['x']
        datay = source['spatial']['y']

        _x = cluster.x-x
        _y = cluster.y-y

        cg_dna = [hash('%1.8f %1.8f'%(g.x,g.y)) for g in cluster.galaxies]
        source_dna = [hash('%1.8f %1.8f'%(datax[i],datay[i])) for i in xrange(len(datax))]

        _in_source = []
        to_move = []


        for i in xrange(len(cg_dna)):
            if cg_dna[i] in source_dna:
                to_move.append(indexOf(source_dna,cg_dna[i]))
                _in_source.append(cluster.galaxies[i])

        if len(to_move) != len(cluster.galaxies):
            print 'Warning - incomplete detection of cluster galaxies in source catalogue!'

        if len(to_move) > 0:
            to_move = numpy.array(to_move)

            if duplicate:
                original = {}
                for key in source.keys():
                    original[key] = []
                    if type(source[key]) == type({}) and key != 'shortcuts':
                        original[key] = {}
                        for sub_key in source[key].keys():
                            if source[key][sub_key] != None:
                                original[key][sub_key] = deepcopy(source[key][sub_key][to_move])
                    else:
                        if source[key] != None and key != 'shortcuts':
                            original[key] = deepcopy(source[key][to_move])


            if displacement == 'spatial':
                source['spatial']['x'][to_move] = source['spatial']['x'][to_move]-_x
                source['spatial']['y'][to_move] = source['spatial']['y'][to_move]-_y
            elif displacement == 'colour':
                if colour == 'auto':
                    col = cluster.parent.selection['defaults']['colour']
                else:
                    col = colour
                mag = cluster.parent.selection['colours'][col]['magnitude']
                #CMR: x=> mag, y=>colour
                source['magnitudes'][mag][to_move] = source['magnitudes'][mag][to_move]+x
                source['colours'][col][to_move] = source['colours'][col][to_move]+y
                if shunt != False:
                    #assume shunt is in arcmins:
                    
                    source['spatial']['x'][to_move] = source['spatial']['x'][to_move]-radians(shunt/60.0)
                    source['spatial']['y'][to_move] = source['spatial']['y'][to_move]-radians(shunt/60.0)
                
            if duplicate:
                for key in original.keys():
                    if type(original[key]) == type({}):
                        for sub_key in original[key].keys():
                            source[key][sub_key] = numpy.array(list(source[key][sub_key])+list(original[key][sub_key]))
                    else:
                        source[key] = numpy.array(list(source[key])+list(original[key]))
                source['x'] = source['spatial']['x']
                source['y'] = source['spatial']['y']

        if tracer:
            if not source.has_key('extras'):
                source['extras'] = {}
            data_range = arange(len(source['spatial']['x']))
            trace = numpy.array([_d in to_move for _d in data_range])
            source['extras']['tracer'] = trace
            


        if do_detection and mode == 'cuts':
            output = parallel_detect(source,det.selection)
            det.addclusters(output)
            return det
        if return_source:
            return source
            
    def hole_test(self,plot=False,rmax=7.0,keep_cuts=True,filling_factor='default'):
        from parallel_vertices import deploy_vertices
        from parallel_qhull import qhull
        from parallel_detect import parallel_detect
        from dataIO import sel_edit,pp_run
        from science import science
        from tools import common_area
        sci = science()
        dependent_functions = (deploy_vertices,qhull,parallel_detect,sel_edit,sci.richness)
        modules = ('science','parallel_detect','clusters','numpy','parallel_qhull','parallel_vertices','vperc','time','detections')
        run_function = 'parallel_detect'
        args = [{'source[i].cuts':None},{'source[i].selection':None}]            
        job_server = None
        ncpus = 16
        
        from tools import objectify,hashify
        from clusters import galaxy,cluster
        from clusters import galaxies_within_r
        nclusters = len(self.clusters)

        try:
            data = self.cuts
            xx = data['x']
        except:
            keep_cuts = False

        if not keep_cuts:
            self.get_input()


        data = self.cuts
        data_gals = objectify(data)

        rmax = radians(rmax/60.0)

        dummy_sel = deepcopy(self.selection)
        det = detections(None,dummy_sel,dummy=True)

        survey_area = self.selection['detection']['detection_area_radians']
        dx = self.selection['limits']['xmax']-self.selection['limits']['xmin']
        dy = self.selection['limits']['ymax']-self.selection['limits']['ymin']

        n = 100.0
        deltax = dx/n
        deltay = dy/n
        gdy = arange(self.selection['limits']['ymin'],self.selection['limits']['ymax']+min(deltay,deltax)/2.0,min(deltay,deltax))
        gdx = arange(self.selection['limits']['xmin'],self.selection['limits']['xmax']+min(deltay,deltax)/2.0,min(deltay,deltax))

        print 'Generating grid galaxies'
        grid_galaxies = []
#        for _x in gdx:
#            for _y in gdy:
#                grid_galaxies.append(galaxy(-1,_x,_y,dummy=True))
        
        print 'Created grid_galaxies'
        grid_in_void = []
        import random
        from math import pi
        dets = []
        det.cuts = data
        dets.append(deepcopy(det))

        if filling_factor == 'default':
            filling_factor = arange(0.01,0.1,0.02)

        all_voids = []
        for ff in filling_factor:
            area = 0.0
            voids = []
            while float(area)/float(survey_area) < ff:
                xpos = self.selection['limits']['xmin'] + random.random()*dx
                ypos = self.selection['limits']['ymin'] + random.random()*dy
                radius = random.random()*rmax
                centroid = galaxy(-1,xpos,ypos,dummy=True)

                voids.append([xpos,ypos,radius])
                total_void_area = sum([pi*(__v[2]**2) for __v in voids])
                overlaps = common_area(voids)

                if (not overlaps < 0) and (not overlaps > 0) and not overlaps == 0:
                    print 'Nan alert'
                    laps = common_area(voids)

                area = total_void_area - overlaps



            eliminated = []
            for v in voids:
                centroid = galaxy(-1,v[0],v[1],dummy=True)
                inside = galaxies_within_r(centroid,data_gals,r=v[2],shape='circle')
                eliminated = eliminated + inside


            all_voids.append(voids)
            eliminated = unique(eliminated)

            if plot:
                self.plot()
                for v in voids:
                    cir = plt.Circle((v[0],v[1]), radius=v[2],fc='gray',alpha=0.5,lw=2.0,ec='gray')
                    plt.gca().add_patch(cir)
                ff0 = deepcopy(float(area)/float(survey_area))
                plt.title('Target ff = %f, actual = %f'%(ff,ff0))

            new_cuts = []
            for g in data_gals:
                if not g in eliminated:
                    new_cuts.append(g)

            cuts = deepcopy(new_cuts)
            cuts = hashify(new_cuts)
            det.cuts = cuts
            dets.append(deepcopy(det))
        if len(dets) > 1:
            completeness = []
#            cls = parallel_detect(dets[1].cuts,dets[1].selection)


            completed,job_server = pp_run(dets,run_function,args,dependent_functions,modules,ncpus=ncpus,idle_delay=360.0,verbose=False,server=job_server)
            for i in xrange(len(dets)):
                dets[i].addclusters(completed[i])

            d0 = dets[0]
            d0_gals = d0.cluster_galaxies

            dets = dets[1:]
            for i in xrange(len(all_voids)):
                #what could have been detected?
                undetectable = []

                hidden_galaxies = []
                for v in all_voids[i]:
                    centroid = galaxy(-1,v[0],v[1],dummy=True)
                    inside = galaxies_within_r(centroid,d0_gals,r=v[2],shape='circle')
                    hidden_galaxies = hidden_galaxies + inside
                parents = [g.parent for g in hidden_galaxies]

                uparents = unique(parents)
                for c in uparents:
                    gals_in_void = float(countOf(parents,c))
                    total_gals = float(len(c.galaxies))
                    if gals_in_void/total_gals >= 0.5:
                        undetectable.append(c)
                                      
                detectable = []
                for c in d0.clusters:
                    if c not in undetectable:
                        detectable.append(c)
                        
                detectable_dna = [hash('%1.8f %1.8f'%(c.bcg.x,c.bcg.y)) for c in detectable]
                
                #what was detected?
                retrieved_dna = [hash('%1.8f %1.8f'%(g.x,g.y)) for g in dets[i].cluster_galaxies]

                match = []
                for dna in detectable_dna:
                    if dna in retrieved_dna:
                        match.append(dna)

                c = float(len(match))/float(len(detectable))
                completeness.append(c)

            return [filling_factor,completeness]

    def determine_overlaps(self,detection,colour='default'):
        from copy import deepcopy
        from detections import detections as det
        if colour == 'default':
            colour = self.selection['defaults']['colour']
            
        try:
            zones = self.zones
            self_clusters = zones['always']
            self_det = det(None,deepcopy(self.selection),dummy=True)
            self_det.addclusters(self_clusters,verbose=False)
        except:
            self_det = self

        try:
            zones = detection.zones
            det_clusters = zones['always']
            other_det = det(None,deepcopy(detection.selection),dummy=True)
            other_det.addclusters(det_clusters,verbose=False)
        except:
            other_det = detection
        

        from science import science
        sci = science()
        alldets = [detection,self]
        xmin = min([d.selection['limits']['xmin'] for d in alldets])
        xmax = max([d.selection['limits']['xmax'] for d in alldets])
        ymin = min([d.selection['limits']['ymin'] for d in alldets])
        ymax = min([d.selection['limits']['ymax'] for d in alldets])

        shared_zone_clusters = []

        zones = {'self':{},'adjacent':{}}
        zones['self']['always'] = []
        zones['self']['choose'] = []
        zones['self']['never'] = []

        zones['adjacent']['always'] = []
        zones['adjacent']['choose'] = []
        zones['adjacent']['never'] = []

        geometry = {}

        d = detection
        d = other_det
        if 1:
            dx = d.selection['limits']['xmax']-d.selection['limits']['xmin']
            dy = d.selection['limits']['ymax']-d.selection['limits']['ymin']
            xy = (d.selection['limits']['xmin'],d.selection['limits']['ymin'])
            self_overlaps,d_overlaps,overlap_region = sci.cherrypick([self,d],colour=colour,two_det=True)
            if len(self_overlaps) > 0:
                shared_zone_clusters = shared_zone_clusters + self_overlaps
            geometry['overlap_region'] = deepcopy(overlap_region)
            dr = overlap_region['width']/3.0
            r1 = {'xy':deepcopy(overlap_region['xy']),'width':overlap_region['width'],'height':overlap_region['height']}
            overlap_region['xy'][0] = overlap_region['xy'][0] + dr
            r2 = {'xy':deepcopy(overlap_region['xy']),'width':overlap_region['width'],'height':overlap_region['height']}
            overlap_region['xy'][0] = overlap_region['xy'][0] + dr
            r3 = {'xy':deepcopy(overlap_region['xy']),'width':overlap_region['width'],'height':overlap_region['height']}


            sub_regions = {'r1':{'geometry':r1},'r2':{'geometry':r2},'r3':{'geometry':r3}}


            print r1
            print r2
            print r3

            if (d.selection['limits']['xmax'] > self.selection['limits']['xmin']) and (d.selection['limits']['xmin'] < self.selection['limits']['xmin']):
                print 'case:   [   d   { ]   self   }'

                sub_regions['r3']['colour'] = 'red'
                sub_regions['r1']['colour'] = 'green'
#                self_filter = sci.filter(self,'x',r3['xy'][0],self.selection['limits']['xmax'])
#                d_filter = sci.filter(d,'x',d.selection['limits']['xmin'],r1['xy'][0]+dr)
#                sz_self = sci.filter(self,'x',r2['xy'][0],r3['xy'][0])
#                sz_d = sci.filter(d,'x',r2['xy'][0],r3['xy'][0])



#                zones['self']['always'] = sci.filter(self_det,'x',r3['xy'][0],self.selection['limits']['xmax']).clusters
                zones['self']['never'] = sci.filter(self_det,'x',d.selection['limits']['xmin'],r2['xy'][0]).clusters
                zones['self']['choose'] = sci.filter(self_det,'x',r2['xy'][0],r3['xy'][0]).clusters

#                zones['adjacent']['always'] = sci.filter(other_det,'x',d.selection['limits']['xmin'],r2['xy'][0]).clusters
                zones['adjacent']['never'] = sci.filter(other_det,'x',r3['xy'][0],self_det.selection['limits']['xmax']).clusters
                zones['adjacent']['choose'] = sci.filter(other_det,'x',r2['xy'][0],r3['xy'][0]).clusters

            else:
                print 'case:   {   self   [ }   d    ]'
                sub_regions['r3']['colour'] = 'green'
                sub_regions['r1']['colour'] = 'red'

#                zones['self']['always'] = sci.filter(self_det,'x',self.selection['limits']['xmin'],r2['xy'][0]).clusters
                zones['self']['never'] = sci.filter(self_det,'x',r3['xy'][0],other_det.selection['limits']['xmin']).clusters
                zones['self']['choose'] = sci.filter(self_det,'x',r2['xy'][0],r3['xy'][0]).clusters

#                zones['adjacent']['always'] = sci.filter(other_det,'x',r3['xy'][0],d.selection['limits']['xmax']).clusters
                zones['adjacent']['never'] = sci.filter(other_det,'x',self_det.selection['limits']['xmin'],r2['xy'][0]).clusters
                zones['adjacent']['choose'] = sci.filter(other_det,'x',r2['xy'][0],r3['xy'][0]).clusters

            excluded = []
            for cc in self_det.clusters:
#                if (not cc in zones['self']['always']) and (not cc in zones['self']['choose']):
                if (not cc in zones['self']['never']) and (not cc in zones['self']['choose']):
                    zones['self']['always'].append(cc)
#                    zones['self']['never'].append(cc)

            excluded = []
            for cc in other_det.clusters:
#                if (not cc in zones['adjacent']['always']) and (not cc in zones['adjacent']['choose']):
                if (not cc in zones['adjacent']['never']) and (not cc in zones['adjacent']['choose']):                    
#                    zones['adjacent']['never'].append(cc)
                    zones['adjacent']['always'].append(cc)
            
            geometry['sub_regions'] = sub_regions


            self.re_posess()
            detection.re_posess()
            return zones,geometry



    def plot_overlaps(self,detection,colour='default',source_colour='red',adjacent_colour='green',overlap_colour='yellow',border=0.02):
        if colour == 'default':
            colour = self.selection['defaults']['colour']


        zones,geometry = self.determine_overlaps(detection,colour=colour)
        overlap_region = geometry['overlap_region']
        sub_regions = geometry['sub_regions']
        
        from science import science
        d = detection
        alldets = [detection,self]
        xmin = min([det.selection['limits']['xmin'] for det in alldets])
        xmax = max([det.selection['limits']['xmax'] for det in alldets])
        ymin = min([det.selection['limits']['ymin'] for det in alldets])
        ymax = min([det.selection['limits']['ymax'] for det in alldets])

        fig = plt.figure()
        sp = plt.subplot(211, axisbg='k',aspect='equal')

        dx = self.selection['limits']['xmax']-self.selection['limits']['xmin']
        dy = self.selection['limits']['ymax']-self.selection['limits']['ymin']
        xy = (self.selection['limits']['xmin'],self.selection['limits']['ymin'])

        # self detection box
        rec = plt.Rectangle(xy,dx,dy,fc=source_colour,alpha=0.4,lw=2.0,ec=source_colour,ls='dashed')
        plt.gca().add_patch(rec)

        dx = d.selection['limits']['xmax']-d.selection['limits']['xmin']
        dy = d.selection['limits']['ymax']-d.selection['limits']['ymin']
        xy = (d.selection['limits']['xmin'],d.selection['limits']['ymin'])

        # adjacent detection box
        rec = plt.Rectangle(xy,dx,dy,fc=adjacent_colour,alpha=0.4,lw=2.0,ec=adjacent_colour,ls='dashed')
        plt.gca().add_patch(rec)

        # overlapping region between the two
        rec = plt.Rectangle(overlap_region['xy'],overlap_region['width'],overlap_region['height'],fc=overlap_colour,alpha=0.4,lw=0.0,ec=overlap_colour,ls='dashed')
        plt.gca().add_patch(rec)

        dr = overlap_region['width']/3.0
        #sub regions - r1:

        if sub_regions['r1']['colour'] == 'red':
            sub_regions['r1']['colour'] = source_colour
            sub_regions['r3']['colour'] = adjacent_colour
        else:
            sub_regions['r1']['colour'] = adjacent_colour
            sub_regions['r3']['colour'] = source_colour
            
            
        rec = plt.Rectangle(sub_regions['r1']['geometry']['xy'],dr,overlap_region['height'],fc=sub_regions['r1']['colour'],alpha=0.4,lw=0.0,ec=sub_regions['r1']['colour'],ls='dashed')
        plt.gca().add_patch(rec)

        #sub regions - r3:
        rec = plt.Rectangle(sub_regions['r3']['geometry']['xy'],dr,overlap_region['height'],fc=sub_regions['r3']['colour'],alpha=0.4,lw=0.0,ec=sub_regions['r3']['colour'],ls='dashed')
        plt.gca().add_patch(rec)

        #now plot the clusters

        plot_clusters(zones['self']['always'],rings=True,centroids=True,passback=True,limits=self.selection,colour='!'+source_colour)
        plot_clusters(zones['self']['choose'],rings=True,centroids=True,passback=True,limits=self.selection,colour='!'+overlap_colour)

        plot_clusters(zones['adjacent']['always'],rings=True,centroids=True,passback=True,limits=self.selection,colour='!'+adjacent_colour)
        plot_clusters(zones['adjacent']['choose'],rings=True,centroids=True,passback=True,limits=self.selection,colour='!'+overlap_colour)
        
        border = dx*border
        plt.xlim(xmin-border,xmax+border)
        plt.ylim(ymin-border,ymax+border)

        sp = plt.subplot(212, axisbg='k',aspect='equal')
        rec = plt.Rectangle((xmin,ymin),xmax-xmin,ymax-ymin,fc='None',alpha=0.0,lw=1.0,ec='gray',ls='dashed')
        plt.gca().add_patch(rec)

        xpoints = [sub_regions['r%d'%(ii)]['geometry']['xy'][0] for ii in range(1,4)]
        # (I rock!)
        xpoints.append(xpoints[-1]+dr)

        for xx in xpoints:
            plt.plot([xx,xx],[ymin,ymax],color='gray',ls='dashed',lw=1.0)

        plot_clusters(self.clusters,centroids=True,rings=True,colour='!'+source_colour,passback=True)
        plot_clusters(d.clusters,centroids=True,rings=True,colour='!'+adjacent_colour,passback=True)
        plt.xlim(xmin-border,xmax+border)
        plt.ylim(ymin-border,ymax+border)
        plt.subplots_adjust(top=0.99,bottom=0.02,hspace=0.04)

    def convert_to_external(self,filename='external.dat',wrap=False):
        from math import pi
        string_buffer=[]
        for c in self.clusters:
            cx,cy = c.x,c.y
            if wrap == False:
                if cx < 0:
                    cx = cx + 2.0*pi
            try:
                string = '%f %f %f\n'%(degrees(cx),degrees(cy),degrees(c.extent*60.0))
            except:
                string = '%f %f %f\n'%(degrees(cx),degrees(cy),degrees(c.max_radius()*60.0))
            string_buffer.append(string)
        if '/' in filename == False:
            filename = './'+filename
        file = open(filename,'w')
        file.writelines(string_buffer)
        file.close()
        print 'External data written to %s'%(filename)
        return
        
      

    def convert(self,verbose=True,all=False):
        selection = {}
        selection['mock_data'] = {}
        selection['colours'] = {}
        selection['defaults'] = {}
        selection['limits'] = {}
        selection['overrides'] = {}
        selection['detection'] = {}
        selection['io'] = {}

        keywords = {}
        keywords['mock_data'] = ['halomass']
        keywords['colours'] = ['width','slope','intercept']
        keywords['defaults'] = ['this_colour']
        keywords['limits'] = ['xmax','xmin','ymax','ymin','zmin','zmax','cmr_range','rlimit']
        keywords['overrides'] = []
        keywords['detection'] = ['dataset','detection_area','detection_area_radians','mean_density','ngal','probthresh',\
                                 'rho_crit','smallest_cluster','boundary','subset_index','metric','empty']
        keywords['io'] = ['dataset','observed_data','wrap','thumb_server','sourcefile']

        source_keys = self.selection.keys()
        added = []

        for key in self.selection.keys():
            for search_key in keywords.keys():
                for word in keywords[search_key]:
                    if key == word:
                        if key == 'this_colour':
                            selection[search_key]['colour'] = self.selection[key]
                            added.append(key)
                        elif key == 'rlimit':
                            selection[search_key]['magnitudes'] = {'r':{'limit':self.selection[key]}}
                            added.append(key)
                        else:
                            selection[search_key][key] = self.selection[key]
                            added.append(key)
                    elif '_' in key and key.split('_')[1] == word:
                        if selection['colours'].has_key(key.split('_')[0]):
                            selection['colours'][key.split('_')[0]][key.split('_')[1]] = self.selection[key]
                            added.append(key)
                        else:
                            selection['colours'][key.split('_')[0]] = {}
                            selection['colours'][key.split('_')[0]][key.split('_')[1]] = self.selection[key]
                            added.append(key)
                    elif type(self.selection[key]) == type({}) and search_key in keywords.keys() and word in self.selection[key].keys():
                        selection[key][word] = self.selection[key][word]
                        added.append(word)
                        

        if 'gr' in selection['colours']:
            selection['detection']['sequence'] = ['gr','ri','iz']
            for key in selection['colours'].keys():
                if key not in selection['detection']['sequence']:
                    del selection['colours'][key]
                    try:
                        del selection['limits']['cmr_range'][key]
                    except:
                        pass
        

        for colour in selection['colours'].keys():
            if not selection['colours'][colour].has_key('fitfuncs'):
                selection['colours'][colour]['fitfuncs'] = {}
            


        if selection['detection'].has_key('sequence') and len(selection['detection']['sequence'][0]) == 2:
            for colour in selection['detection']['sequence']:
                selection['colours'][colour]['magnitude'] = colour[-1]

        missing = []
        for key in source_keys:
            if (not key in added) and (not key in selection.keys()):
                missing.append(key)
        if verbose:
            print 'Items not added to new selection : %s'%(str(missing))

        self.selection = selection

        if all:
            #1 check neighbour galaxies
            if self.field_galaxies != None:
                neighbours = []
                for g in self.field_galaxies:
                    _g = g.convert()
                    if _g:
                        neighbours.append(_g)
                    else:
                        neighbours.append(g)
                self.field_galaxies = neighbours
                        
                
            #2 check clusters
            clusters = []
            for c in self.clusters:
                _c = c.convert()
                if _c:
                    clusters.append(_c)
                else:
                    clusters.append(c)
            if len(clusters) > 0:
                self.addclusters(clusters,verbose=False,perfect=False)
            else:
                self.clusters = []


    def get_input(self,data=None,perfect=False,pick_list=None):
        from dataIO import make_mask,getData,slice_data
        from clusters import photometric_selection
        from detections import detections as det
        import numpy
        import time
#        pick_list = ['x','y','r']
#        pick_list = pick_list + self.selection['detection']['sequence']

        try:
            observed = self.selection['io']['observed_data']
        except:
            observed = True

        if (data == None or data == False):
            try:
                data = self.source_data
                xx = data['spatial']['x']
            except:            
                wrap = False
                if self.selection['limits']['xmin'] > 2.0 and 's82' in self.selection['io']['dataset'] and 'blk' not in self.selection['io']['dataset']:
                    wrap = True
                    self.wrap()
                if self.selection['limits']['xmin'] < 0:
                    wrap = True
                data = getData(self.selection['io']['sourcefile'],obs=observed,dataset=self.selection['io']['dataset'],selection=self.selection,wrap=wrap)


        if self.selection.has_key('overrides'):
            if self.selection['overrides'].has_key('project_cluster'):
                self.source_data = data
                pc = self.selection['overrides']['project_cluster']
                _cluster = pc['cluster']
                _x = pc['x']
                _y = pc['y']
                print 'Projecting cluster... @ ra=%f,dec=%f'%(_x,_y)
                data = self.project_to(_cluster,_x,_y,do_detection=False,return_source=True,keep_cuts=True,mode='data')
                self.source_data = data
                
            if self.selection['overrides'].has_key('implant_cluster'):
                self.source_data = data
                ic = self.selection['overrides']['implant_cluster']
                _cluster = ic['cluster']
                _x = ic['x']
                _y = ic['y']
                print 'Implanting cluster... @ ra=%f,dec=%f'%(_x,_y)
                data = self.implant_cluster(_cluster,_x,_y,do_detection=False,return_source=True,keep_cuts=True,mode='data')
                self.source_data = data

#        del self.selection['grid']
#        print self.selection['colours']
#        print self.selection['limits']
        
#        print self.selection['colours']
#        raise SystemExit(0)
        mask = photometric_selection(data,self.selection)
        arr = numpy.array([])
        new = {}
        for key in data.keys():
            new[key] = {}
            if key == 'shortcuts':
                continue
            elif data.has_key('shortcuts'):
                if key in data['shortcuts']:
                    continue
            for sub_key in data[key].keys():
                if type(data[key][sub_key]) != type(None) and type(data[key][sub_key]) == type(arr):
                    if pick_list == None:
                        new[key][sub_key] = numpy.take(data[key][sub_key],mask)
                    elif sub_key in pick_list:
                        new[key][sub_key] = numpy.take(data[key][sub_key],mask)
                else:
                    new[key][sub_key] = None


        if new.has_key('shortcuts'):
            for key in data['shortcuts']:
                new['shortcuts'][key] = {}
                for _key in new.keys():
                    if _key != 'shortcuts' and _key not in new['shortcuts'].keys():
                        for sub_key in new[_key].keys():
                            if sub_key == key:
                                new[sub_key] = new[_key][sub_key]
        
            

        self.cuts = new
        string = '\t#Made c('
        colours = ''
        ints = '('
        for col in self.selection['detection']['sequence']:
            colours = colours + '%s,'%(col)
            ints = ints + '%7.5f,'%(self.selection['colours'][col]['intercept'])
        colours = colours[:-1]+')'
        ints = ints[:-1]+')'
        
        dummy = det(self.cuts,self.selection)

        self.selection = deepcopy(dummy.selection)
        ngal = len(self.cuts['spatial']['x'])

        if 'mock' in self.selection['io']['dataset'] and perfect:
            print 'Getting perfect clusters'
            outlier_metric = self.selection['mock_data']['outlier_metric']
            halomass = self.selection['mock_data']['halomass']
            
            self.perfect(data,halomass=halomass,outlier_metric=outlier_metric)

        string = '%s%s=%s detection class with %d galaxies'%(string,colours,ints,ngal)
        print string

    def convert_from_external(self,filename='',wrap=True):
        from clusters import cluster,galaxy
        base = self
        base.compare(filename,name='internal',wrap=wrap)
        
        cluster_list = []
        data = base.external_data['internal']
        x,y = data['x'],data['y']
        if data.has_key('r'):
            r = data['r']
        else:
            r = [0.0 for i in x]

        for i in xrange(len(x)):
            _gal = galaxy(-1,x[i],y[i],dummy=True)
            _cl = cluster([_gal],0,dummy=True)
            _cl.extent = radians(r[i]/60.0)
            cluster_list.append(_cl)
        del base.external_data['internal']
        return cluster_list

    def boundary_check(self,ra,dec,constraint=None,status=None,pp=False,ncpu=16):
        """\
        det.boundary_check(self,ra,dec,constraint=None,status=None,pp=False) => []
        \t determines if dataset is inside this darb boundary
        
        sampling_freq   :  
        ra              : a *list* of ra to check
        dec             : corresponding *list* of dec
        constraint      : list of boolean elements to pre-filter positions
        status          : progress update
        pp              : run parallel
        """

        from clusters import point_inside_cell
        import numpy

        if pp:
            status = None

        if not self.selection.has_key('boundary'):
            print 'No boundary data - returning'
            return

        else:
            if constraint == None:
                constraint = numpy.ones(len(ra),dtype='bool')
            islands = []
            holes = []
            if self.selection.has_key('grid') and self.selection['grid']['self'].has_key('boundary_class') and sum([len(v) for v in self.selection['grid']['self']['boundary_class'].values()]) > 0:
                print 'Using boundary_class shortcuts - worry here about self.selection[boundary] vs self.selection[grid][self][boundary]'
                for j in self.selection['grid']['self']['boundary_class']['island']:
#                    islands.append(self.selection['grid']['self']['boundary'][j])
                    islands.append(self.selection['boundary'][j])
                for j in self.selection['grid']['self']['boundary_class']['hole']:
#                    holes.append(self.selection['grid']['self']['boundary'][j])
                    holes.append(self.selection['boundary'][j])
                    
            else:
                from clusters import point_inside_cell
                var_type = __builtins__['type']          #### horrid! let that be a lesson to me - *dont* override built-ins from namespace.
                partitions = []
                if var_type(self.selection['boundary'][0][0]) == var_type([]):
                    partitions = self.selection['boundary']
                    islands = []
                else:
                    islands = [self.selection['boundary']]
                holes = []

                if len(partitions) > 0:
#                    from copy import deepcopy
#                    plen = [len(p) for p in partitions]
#                    splen = deepcopy(plen)
#                    splen.sort()
#                    partitions = [partitions[indexOf(plen,s)] for s in plen]
                    partitions_rev = partitions

                    if status:
                        from dataIO import Progress
                        status = Progress(len(partitions))

                    for p in partitions:
                        if status: 
                            status = Progress(len(partitions_rev))
                        for p2 in partitions_rev:
                            if p != p2 and (p2 not in holes):
                                is_hole = False
                                points = iter(p)
                                hole_status = []
                                while True:
#                                while (is_hole == False):
                                    try:
                                        xy = points.next()
                                        is_hole = point_inside_cell(xy[0],xy[1],p2)
                                        hole_status.append(is_hole)
                                    except StopIteration:
                                        break


                                if not (False in hole_status):
                                    holes.append(p)
                                


                                junk = """\
                                #this can be sped up with a while + break on 1st True....
                                test = [point_inside_cell(xy[0],xy[1],p2) for xy in p]
                                if countOf(test,True) == len(test):
                                    # this means p is inside p2    => a hole
                                    holes.append(p)
                                elif countOf(test,False) == len(test):
                                    # this means p is not inside p2   => still perhaps a hole though
                                    pass
                                else:
                                    #this suggests p straddles p2 - hmm?!
                                    print 'Panic!!'
                    """

                    for p in partitions:
                        if not p in holes:
                        #must be an island
                            islands.append(p)
                            
                if self.selection.has_key('grid') and var_type(self.selection['boundary'][0][0]) == var_type([]) and len(islands) > 1:
                    print 'setting boundary_class shortcuts for future'
                    self.selection['grid']['self']['boundary_class'] = {}
                    self.selection['grid']['self']['boundary_class']['island'] = []
                    self.selection['grid']['self']['boundary_class']['hole'] = []
                    for region in islands:
                        j = indexOf(self.selection['boundary'],region)
                        self.selection['grid']['self']['boundary_class']['island'].append(j)
                    for region in holes:
                        j = indexOf(self.selection['boundary'],region)
                        self.selection['grid']['self']['boundary_class']['hole'].append(j)
                        
            if status:
                from dataIO import Progress
                print 'Checking data locations... \n'
                status = Progress(len(ra))

            in_survey = numpy.zeros(len(ra),dtype='bool')
            ncpus = ncpu
            if pp:
                n = int(len(ra)/float(ncpus))+1
                
                master_count = iter(arange(len(ra)))
                count = 0
                ras = []
                decs = []
                _ra = []
                _dec = []
                while True:
                    try:
                        i = master_count.next()
                        if count == n:
                            decs.append(_dec)
                            ras.append(_ra)
                            count = 0
                            _ra = []
                            _dec = []
                        _ra.append(ra[i])
                        _dec.append(dec[i])
                        count += 1
                    except StopIteration:
                        decs.append(_dec)
                        ras.append(_ra)
                        break
                print 'Divided... now conquer'
                args = [{'source[i]':None},{'ra[i]':ras},{'dec[i]':decs}]
                ncpus = 16
                dependent_functions = ()
                modules = ('science','parallel_detect','clusters','numpy','parallel_qhull','parallel_vertices','vperc','time','detections')
                run_function = 'boundary_check'
                job_server = None

                from copy import deepcopy
                from dataIO import pp_run
            
                core = self.clone()
                dets = [deepcopy(core) for d in decs]

                completed,job_server = pp_run(list(dets),run_function,args,dependent_functions,modules,ncpus=ncpus,idle_delay=360.0,verbose=True,server=job_server,progress=False)
                print 'PP job completed, merging'
                output = []
                for c in completed:
                    output = output + list(c)
                output = numpy.array(output)
                output = (output & constraint)
                return output

            if len(islands) == 0 and len(holes) > 0:
                #this is assumed a field covering the entire range, with holes inside => reverse case
                in_survey = numpy.ones(len(ra),dtype='bool')
                
            
            for i in xrange(len(ra)):
                _x,_y = ra[i],dec[i]
                inside = False
                for island in islands:
                    inside = point_inside_cell(_x,_y,island)
                    if inside == True:
                        break
                in_hole = False
                if inside == True:
                    for hole in holes:
                        in_hole = point_inside_cell(_x,_y,hole)
                        if in_hole == True:
                            break
                if (in_hole == False) and (inside == True):
                    in_survey[i] = not in_survey[i]

                if status: 
                    status.update(0)


            junk = """\
            for i in xrange(len(ra)):
                if constraint[i]:
                    _x,_y = ra[i],dec[i]
                    if len(holes) > 0:
                        inside = [point_inside_cell(_x,_y,h) for h in holes]
                        if True in inside:
                            in_survey[i] = False
                            continue
                    if len(islands) > 0:
                        inside = [point_inside_cell(_x,_y,j) for j in islands]
                        if True in inside:
                            in_survey[i] = True
                if status: 
                    status.update(0)

            """ 
        return in_survey
        
    def ingest_weightmap(self,map):
        all_galaxies = self.cluster_galaxies+self.associated_member_galaxies
        all_galaxies = self.cluster_galaxies
        try:
            all_galaxies = all_galaxies+self.associated_member_galaxies
        except:
            pass
        for g in all_galaxies:
            g.wsel = map.get_weight(degrees(g.x),degrees(g.y))


    def compare(self,catalogue,plot=False,passback=False,colour='auto',name='unspecified',clip=True,wrap=False,exact=True,rings=True,points=True,shape='auto',type='clusters',style='solid',import_all=False,default_scale=8.0,check_border=True,lw=1.0,radial=True):
        """\
        detections.compare(catalogue,plot=False,passback=False,colour='auto',name='unspecified',clip=True,wrap=False,exact=True,rings=True,points=True,shape='auto',type='clusters',style='solid',import_all=False,default_scale=8.0,check_border=True,lw=1.0,radial=True) -> void
    
        \t Import a set of alternative detections with which to compare this detection with. These can
        be clusters themselves, or galaxies

        catalogue     : Can be either a string (pointing to an ASCII file, a .fits file), or a list of positions
                        Formats
                        ASCII: * lists of    RA DEC [z or r]
                               * IAU-style positions (e.g. J035020.6+133741.4)
                        FITS : * file must at least contain columns with 'ra' and 'dec'. 
                                 Will import 'z' and 'Name' if they exist
                        list : * supply as [[ra0,dec0,|z0|],[ra1,dec1,|z1|]....]
        plot          : make a plot of the data when ingested?
        passback      : allow further additions to this figure after plotting?
        colour        : colour of the symbols plotted for the data
        name          : an identifier for the data ingested (e.g. Abell, BCS)
        clip          : limit input dataset to the boundary (approx) of the detection class?
        wrap          : deal with wrap-around for Stripe 82 (where the RA is -ve in some treatments for data contiguity)
        exact         : set the limits (see clip) to precisely those of the detection class
        rings         : plot data as circles (if plot=True)
        points        : plot a marker at the detection location
        type          : what is the data being ingested? ['clusters'|'galaxy']
        style         : what line-style should the line plots be for this dataset?
        import_all    : override all constraints (postitional etc) and import everything supplied
        default_scale : if no size data (e.g. distance) provided, what is the dummy angular size?
        check_border  : if the detection class has been processed with DARB, check entried lie within data boundary
        lw            : plotting line-width (if plot=True)
        radial        : check data for radial (ie distance) data (this is normally automated anyway)


        """



        from dataIO import load
        import numpy
        from math import atan,degrees
        redshift = False
        zz = None
        cl_names = None
        var_type = __builtins__['type']
        extras = None
        e_names = []
        if (var_type(catalogue) == var_type('str')) and ('.fit' in catalogue):
            non_extras = ['ra','dec','Name','z']
            import fits_table as fits
            tbl = fits.fits_table(catalogue)
            ra = tbl.get_column('ra')
            dec = tbl.get_column('dec')
            if 'z' in tbl.names:
                redshift = True
                zz = tbl.get_column('z')
                rad = tbl.get_column('z')

            if 'Name' in tbl.names:
                cl_names = numpy.array(tbl.get_column('Name'))
                print 'Got names'

            for _name in tbl.names:
                if not _name in non_extras:
                    if not extras:
                        extras = {}
                    e_names.append(_name)
                    extras[_name] = tbl.get_column(_name)
            if len(e_names) > 0:
                e_names.sort()

# in catalogue and var_type(load(catalogue,verbose=False)) == var_type(self):



        elif (var_type(catalogue) == var_type(self)) or ((var_type(catalogue) == var_type('str')) and ('.dat' in catalogue) and var_type(load(catalogue,verbose=False)) == var_type(self)):

            if var_type(catalogue) == var_type('str'):
                catalogue = load(catalogue)
            #this is an ORCA detection file?
            print 'ORCA det!'
            from math import degrees
            ra = numpy.array([degrees(cl.x) for cl in catalogue.clusters])
            dec = numpy.array([degrees(cl.y) for cl in catalogue.clusters])
#            ra = numpy.array([cl.x for cl in catalogue.clusters])
#            dec = numpy.array([cl.y for cl in catalogue.clusters])
#

            rad = numpy.array([degrees(cl.extent)*60.0 for cl in catalogue.clusters])
            zz = numpy.array([cl.z for cl in catalogue.clusters])
            refshift = False

#            return

#            if type(det_compare

        elif var_type(catalogue) != var_type('string') and var_type(catalogue) in [var_type([]),var_type(numpy.array([]))]:
            #inline data as arg!
            rad = None
            e1 = catalogue[0]
            if len(e1) not in [2,3]:
                print 'Inline submission of data must be a list of type [[ra0,dec0,|z0|],[ra1,dec1,|z1|]....]'
                return
            ra = [c[0] for c in catalogue]
            dec = [c[1] for c in catalogue]
            if len(c) == 3:
                rad = numpy.array([c[2] for c in catalogue])




        else:


            try:
                data = numpy.loadtxt(catalogue).transpose()
            except:
                data = numpy.loadtxt(catalogue,dtype='string')
                if 'J' in data[0]:
                    radata = []
                    decdata = []
                    from tools import iau_radec
                    for d in data:
                        radec = iau_radec(d)
                        radata.append(radec[0])
                        decdata.append(radec[1])
                    data = [radata,decdata]
    #                print 'Ingested IAU positions'
    #                return data
                else:
                    print 'File %s is not something I can read at the moment'
                    print 'I accept lists of:  RA DEC [z || r]'
                    print 'RA DEC [z || r]          or'
                    print 'IAU-style positions (e.g. J035020.6+133741.4)'
                    return
            ra = data[0]
            dec =data[1]

            try:
                _gg = len(ra)
            except:
                ra = numpy.array([ra])
                dec = numpy.array([dec])

            rad = None
            try:
                rad = data[2]
                try:
                    _gg = len(rad)
                except:
                    rad = numpy.array([rad])

                if max(data[2]) > 3.0:
                    redshift = False
                elif max(data[2]) < 3.0 and min(data[2]) > 0.0:
                    redshift = True

                elif min(rad) < 0:
                    rad = None
                    

            except:
                rad = None
        
        ra = numpy.array(radians(ra))
        dec = numpy.array(radians(dec))
        offset = numpy.zeros(len(ra),dtype='bool')


        if ('s82' in self.selection['io']['dataset'] and wrap) or (self.selection['limits']['xmin'] < 0 and wrap):
            xform = (ra > 2.0)
            ra[xform] = ra[xform] - 2.0*pi
            offset[xform] = True


        elif self.selection['io'].has_key('wrap') and self.selection['io'].has_key('wrap_threshold'):
            print 'Stripe-type thresholding... converting'
            xform = (ra < self.selection['io']['wrap_threshold'])
            ra[xform] = ra[xform] + 2.0*pi


        ra_min = self.selection['limits']['xmin']-(self.selection['limits']['xmax']-self.selection['limits']['xmin'])
        ra_max = self.selection['limits']['xmax']+(self.selection['limits']['xmax']-self.selection['limits']['xmin'])
        dec_min = self.selection['limits']['ymin']-(self.selection['limits']['ymax']-self.selection['limits']['ymin'])
        dec_max = self.selection['limits']['ymax']+(self.selection['limits']['ymax']-self.selection['limits']['ymin'])

        if exact:
            ra_min = self.selection['limits']['xmin']
            ra_max = self.selection['limits']['xmax']
            dec_min = self.selection['limits']['ymin']
            dec_max = self.selection['limits']['ymax']
            

        if clip:
            constraint = (ra > ra_min) & (ra < ra_max) & (dec > dec_min) & (dec < dec_max)
        else:
            constraint = (dec > dec_min) & (dec < dec_max)


            
        if import_all:
            constraint  = numpy.ones(len(ra),dtype='bool')
            
        if len(constraint) > 1000:
            print 'Dataset %s is quite large - pre-filtering'%(name)
            ra = ra[constraint]
            dec = dec[constraint]
            offset = offset[constraint]
            if zz != None:
                zz = zz[constraint]
            if len(e_names) > 0:
                for _name in e_names:
                    extras[_name] = extras[_name][constraint]

            constraint  = numpy.ones(len(ra),dtype='bool')
            print 'Sample size now %d'%(len(constraint))


        if check_border and self.selection.has_key('boundary'):
            in_survey = self.boundary_check(ra,dec,constraint=constraint)
            constraint = constraint & in_survey
        

        ra = list(ra[constraint])
        dec = list(dec[constraint])
        offset = list(offset[constraint])
        if zz != None:
            zz = list(zz[constraint])
        if len(e_names) > 0:
            for _name in e_names:
                extras[_name] = list(extras[_name][constraint])
            
                


#        return cl_names
        if var_type(cl_names) != var_type(None):
            cl_names = list(cl_names[constraint])

        if self.selection['io'].has_key('wrap') and self.selection['io'].has_key('wrap_threshold') and len(offset) > 0:
            _xform = (numpy.array(ra) > 2.0*pi)
            offset = numpy.array(offset)
            offset[_xform] = True
            

        radials = []
        try:
            _r = rad
        except:
            rad = None


        if rad != None:
            rad = list(rad[constraint])
            from tools import comoving
            #assume this is redshift...

            for r in rad:
                if redshift:
                    size = 0.
                    if r > 0:
                        comov = comoving(r)
                        if comov > 0:
                            tan_size = 1.0/comov
                            size = degrees(atan(tan_size))
                            size = size*60.0

                elif radial:
                    size = r
                else:
                    size = default_scale
                radials.append(size)

        else:
            radials = rad

        if len(radials) == 0:
            radials = numpy.ones(len(dec))
            radials = radials*default_scale
            rad = radials

        if type == 'clusters' or type == 'cluster':
            from clusters import cluster
            data = []
            for i in xrange(len(ra)):
                cl = cluster([],-1,dummy=True)
                cl.x = ra[i]
                cl.y = dec[i]
                cl.wrap_alert = offset[i]
                cl.parent = self
                try:
                    cl.get_iau_name(prefix=name)
                except:
                    pass
                if var_type(cl_names) != var_type(None):
                    cl.iau_name = cl_names[i]
                cl.dna = i
                cl.position = i
                if len(radials) > 0:
                    cl.extent = radians(radials[i]/60.0)
                    cl.z = rad[i]
                    if zz != None:
                        cl.z = zz[i]
                if len(e_names) > 0:
                    cl.extras = {}
                    for _name in e_names:
                        cl.extras[_name] = extras[_name][i]
                        
                cl.galaxies = []
                cl.ngalaxies = 1
                cl.luminosity = 0.0
                cl.probability_density = 100.0
                cl.flux = 100.0
                cl.flux_density,cl.flux_density,cl.survey_mean_flux_density = 10.0,10.0,10.0
                data.append(cl)

        elif 'gal' in type:
            from clusters import galaxy
            data = []
            for i in xrange(len(ra)):
                gal = galaxy(0,ra[i],dec[i],dummy=True)
                gal.lum = 0
                gal.parent = self
                if (zz != None) and len(zz) > 0:
                    gal.specz = zz[i]
                if var_type(cl_names) != var_type(None):
                    gal.name = cl_names[i]
                data.append(gal)
                
        try:
            ed = self.external_data
        except:
            self.external_data = None


        if self.external_data != None:
            i = 2
            while self.external_data.has_key(name):
                name = name + '_%s'%(str(i))
            if len(radials) > 0:
                self.external_data[name] = {'x':ra,'y':dec,'r':radials,'colour':'red','shape':'circle','points':points,'type':type,'data':data,'style':style}
            else:
                self.external_data[name] = {'x':ra,'y':dec,'colour':'red','shape':'circle','points':points,'type':type,'data':data,'style':style}
                radials = numpy.ones(len(dec))
                radials = radials*5.0
                self.external_data[name]['r'] = radials
        else:
            if len(radials) > 0:
                self.external_data = {name:{'x':ra,'y':dec,'r':radials,'colour':'red','shape':'circle','points':points,'type':type,'data':data,'style':style}}
            else:
                self.external_data = {name:{'x':ra,'y':dec,'colour':'red','shape':'circle','points':points,'type':type,'data':data,'style':style}}
                radials = numpy.ones(len(dec))
                radials = radials*5.0
                self.external_data[name]['r'] = radials

                
        self.external_data[name]['lw'] = lw
        if colour != 'auto':
            self.external_data[name]['colour'] = colour
        if shape != 'auto':    
            self.external_data[name]['shape'] = shape
        if rad != None:
            self.external_data[name]['redshift'] = rad

        if plot:
            if colour == 'auto':
                colour = 'r'
            self.extdata_plot(passback=passback,override={'name':self.external_data[name]},colour=colour,rings=rings,points=points)

#        return

    def get_density(self):
        ngals = self.selection['ngal']
        da = self.selection['detection_area']
        print ngals,da
        
        return float(ngals)/float(da)


    def extdata_convert(self,colour='auto',shape='circle',label='all',points=True):
        if colour == 'auto':
            paint_gun = ['r','g','b','y']
            plot_colour = iter(paint_gun)
        else:
            paint_gun = [colour]
            plot_colour = iter(paint_gun)

        ed = {}

        keys = self.external_data.keys()
        if label != 'all' and label in self.external_data.keys():
            keys = [label]
        for key in keys:
            if type(self.external_data[key]) == type({}):
                ed[key] = self.external_data[key]
                continue
            ed[key] = {'x':self.external_data[key][0],'y':self.external_data[key][1],'shape':shape,'points':points}
            if len(self.external_data[key]) > 2:
                    if type(self.external_data[key][-1]) == type(' '):
                        ed[key]['colour'] = self.external_data[key][-1]
                        if len(self.external_data[key]) > 3:
                            ed[key]['r'] = self.external_data[key][-2]
                            
                    else:
                        ed[key]['r'] = self.external_data[key][-1]
                        try:
                            ed[key]['colour'] = plot_colour.next()
                        except StopIteration:
                            if colour == 'auto':
                                paint_gun = ['r','g','b','y']
                                plot_colour = iter(paint_gun)
                            else:
                                paint_gun = [colour]
                                plot_colour = iter(paint_gun)
                                ed[key]['colour'] = plot_colour.next()

            else:
                try:
                    ed[key]['colour'] = plot_colour.next()
                except StopIteration:
                    if colour == 'auto':
                        paint_gun = ['r','g','b','y']
                        plot_colour = iter(paint_gun)
                    else:
                        paint_gun = [colour]
                        plot_colour = iter(paint_gun)
                        ed[key]['colour'] = plot_colour.next()

        if len(self.external_data.keys()) == 1 and key == self.external_data.keys()[0] or label=='all':
            self.external_data = ed
            return ed
        else:
            return ed 

    def extdata_plot(self,passback=True,override=None,colour='auto',lw=0.5,rings=True,points=True):
        _lw = lw

        if colour == 'auto':
            paint_gun = ['r','g','b','y']
            plot_colour = iter(paint_gun)
        else:
            paint_gun = [colour,colour,colour,colour,colour,colour,colour,colour]
            plot_colour = iter(paint_gun)

        if override != None:
            datasource = override
        else:
            datasource = self.external_data

        if datasource == None:
            return
        if passback == False:
            fig = plt.figure()
#            sp = plt.subplot(111, axisbg='black',aspect='equal')
            sp = plt.subplot(111, axisbg='black')
            plt.xlabel('RA')
            plt.ylabel('DEC')

        convert = False
        for d in datasource.keys():
            if type(datasource[d]) == type([]):
                convert = True
        if convert:
            self.extdata_convert()
            print 'Converted... now quitting'
            datasource = self.external_data
        styles = iter(['o','s','D','^'])
        added = []
        for data in datasource:
            ra = datasource[data]['x']
            dec = datasource[data]['y']
            have_radii = datasource[data].has_key('r')
            col = datasource[data]['colour']
            if datasource[data].has_key('lw'):
                _lw = datasource[data]['lw']
            else:
                _lw = lw

            if have_radii:
                radii = datasource[data]['r']
                

            if datasource[data] in added:
                continue

            if datasource[data].has_key('data'):
                style = 'solid'
                if datasource[data].has_key('style'):
                    style = datasource[data]['style']

                if 'gal' in datasource[data]['type']:
                    from detections import plot_galaxies
                    plot_galaxies(datasource[data]['data'],ms='x',colour=datasource[data]['colour'],vertices=False,passback=True)
                        
                else:
                    for c in datasource[data]['data']:
                        c.plot(centroid=(datasource[data]['points'] and points),marker_style='o',colour='!'+datasource[data]['colour'],rings=True,shape=datasource[data]['shape'],ring_ls=style,lw=_lw)

            else:
                print 'Old style'
                if datasource[data]['points'] and points:
                    plt.plot(ra,dec,color=col,marker=styles.next(),ls='None')
                shape = datasource[data]['shape']
                if have_radii:
                    for i in xrange(len(radii)):
                        rad = radians(radii[i]/60.0)
                        if shape == 'circle':
                            cir = plt.Circle((ra[i],dec[i]), radius=rad,fc='None',alpha=1.0,lw=lw,ec=col,ls=style)
                        elif shape == 'square':
                            cir = plt.Rectangle((ra[i]-rad,dec[i]-rad),rad*2.0,rad*2.0,fc='None',alpha=1.0,lw=lw,ec=col,ls=style)
                        plt.gca().add_patch(cir)

        if passback == False:
#            plt.xlim(self.selection['limits']['xmin'],self.selection['limits']['xmax'])
            plt.ylim(self.selection['limits']['ymin'],self.selection['limits']['ymax'])
        
    def boss_ingest(self,filename,update=True,write=False,output='boss_ingested.fits',version=1,compare=False,plot=True,dz=0.05):
        import numpy
        import fits_table as fits
        tbl = fits.fits_table(filename)
        clusters = {}
        for cl in self.clusters:
            clusters[cl.dna] = cl
            cl.bossZ = 0
        all_gals = self.cluster_galaxies+self.associated_cluster_galaxies
        galaxies = {}
        objID_list = []
        for g in all_gals:
            galaxies[g.objID] = g
            objID_list.append(g.objID)

        no_dupes = True
        if len(objID_list) != len(unique(objID_list)):
            no_dupes = False

        objID_list = numpy.array(objID_list)
        
        extra_keys = g.extras.keys()

        test = [True]
        key = 'boss'
        version = 1
        while True in test:
            test = [key in k for k in extra_keys]
            if True in test:
                version += 1
                key = 'boss_v%d'%(version)

        objID = tbl.get_column('gal_objID')
        bossZ = tbl.get_column('Z')
        boss_clid = tbl.get_column('clusterID')


        counter = 0
        for g in all_gals:
            g.extras[key] = numpy.nan

        for id,bz in zip(objID,bossZ):
            if id in galaxies.keys():
                for j in numpy.where(objID_list == id)[0]:
                    all_gals[j].extras[key] = bz
                    all_gals[j].parent.bossZ = 1
                    counter += 1

        print 'Added %d boss redshifts into cluster galaxy list [used key: %s]'%(counter,key)
        if compare:
            old_z = numpy.array([cl.z for cl in self.clusters])
            clusterZ = []
            for bz,b_clid in zip(bossZ,boss_clid):
                clusterZ.append(clusters[b_clid].z)
            clusterZ = numpy.array(clusterZ)
            


        if update:
            old_z = numpy.array([cl.z for cl in self.clusters])
            if plot:
                fig = plt.figure()
                
                bins = arange(0,max(old_z)+0.15,dz)
                junk = self.zdist(plot=True,label='clusterZ',colour='r',bins=bins)

            for cl in self.clusters:
                cl.redshift()

            if plot:
                junk = self.zdist(plot=True,label='clusterZ_boss',colour='b',bins=bins)
                plt.legend(loc=1)
                
            new_z = numpy.array([cl.z for cl in self.clusters])
            dz = old_z-new_z
            fig = plt.figure()
            plt.plot(new_z,dz,'ko')
            plt.ylabel('clusterZ-clusterZ_boss')
            plt.xlabel('clusterZ_boss')

            if write:
                from tools import cat
                out = cat(self,output=output,extras=['bossZ'])
            

            return old_z,new_z
    
    def hyperz_ingest(self,filename,update=True,write=False,output='hyperz_ingested.fits',version=1,compare=False,plot=True,dz=0.05,clear=True,backup=False):
        import numpy
        import fits_table as fits
        from dataIO import Progress
        tbl = fits.fits_table(filename)
        clusters = {}
        for cl in self.clusters:
            clusters[cl.dna] = cl
            cl.hyperZ = 0
        all_gals = self.cluster_galaxies+self.associated_cluster_galaxies
        galaxies = {}
        objID_list = []
        for g in all_gals:
            galaxies[g.objID] = g
            objID_list.append(g.objID)

        no_dupes = True
        if len(objID_list) != len(unique(objID_list)):
            no_dupes = False

        objID_list = numpy.array(objID_list)
        
        extra_keys = g.extras.keys()



        test = [True]
        key = 'hyperz'
        if clear and key in extra_keys:
            if backup == False:
                print 'Removing existing %s key'%(key)
            else:
                print 'Backing up existing key %s to %s_backup'%(key,key)
            for g in all_gals:
                if backup:
                    g.extras['%s_backup'%(key)] = g.extras[key]
                del g.extras[key]
            if compare:
                print 'Determining new clusterZ based on removed data...'
                status = Progress(len(self.clusters))
                for cl in self.clusters:
                    cl.redshift()
                    status.update(0)
                print '...done'
        extra_keys = g.extras.keys()
        version = 1
        while True in test:
            test = [key in k for k in extra_keys]
            if True in test:
                version += 1
                key = 'hyperz_v%d'%(version)
        key = 'hyperz'
        objID = tbl.get_column('objID')
        hyperZ = tbl.get_column('hyperz')
        hyperZ_clid = tbl.get_column('clusterID')


        counter = 0
        for g in all_gals:
            g.extras[key] = numpy.nan

        for id,bz in zip(objID,hyperZ):
            if id in galaxies.keys():
                for j in numpy.where(objID_list == id)[0]:
                    all_gals[j].extras[key] = bz
                    all_gals[j].parent.bossZ = 1
                    counter += 1

        print 'Added %d hyperz redshifts into cluster galaxy list [used key: %s]'%(counter,key)
        if compare:
            missing_clusters = []
            old_z = numpy.array([cl.z for cl in self.clusters])
            clusterZ = []
            for bz,b_clid in zip(hyperZ,hyperZ_clid):
                if clusters.has_key(b_clid):
                    clusterZ.append(clusters[b_clid].z)
                else:
                    missing_clusters.append(b_clid)
            clusterZ = numpy.array(clusterZ)
            if len(missing_clusters) > 0:
                missing_clusters = unique(missing_clusters)
                print 'The following clusters were not in the source data:'
                mcl = ''
                for cl in missing_clusters:
                    mcl = mcl + '%s '%(cl)
                print str(mcl)


        if update:
            old_z = numpy.array([cl.z for cl in self.clusters])
            if plot:
                fig = plt.figure()
                
                bins = arange(0,max(old_z)+0.15,dz)
                junk = self.zdist(plot=True,label='clusterZ',colour='r',bins=bins)

            for cl in self.clusters:
                cl.redshift()

            if plot:
                junk = self.zdist(plot=True,label='clusterZ_hyperz',colour='b',bins=bins)
                plt.legend(loc=1)
                
            new_z = numpy.array([cl.z for cl in self.clusters])
            dz = old_z-new_z
            fig = plt.figure()
            plt.plot(new_z,dz,'ko')
            plt.ylabel('clusterZ-clusterZ_hyperz')
            plt.xlabel('clusterZ_hyperz')

            if write:
                from tools import cat
                extras = []
                if backup:
                    extras = []
                out = cat(self,output=output,extras=extras)
            

            return old_z,new_z
    





    def boss_simple(self,boss_file,dr7_file,maglim=21.1,brightlim=17.0,uber=False,band='i',target=3,extras=None,no_list=False,write=True,supplemental=False,bcg_only=False):
        import numpy
        from detections import detections as det
        bcgs = [c.bcg for c in self.clusters]

        dummy = det(None,deepcopy(self.selection),dummy=True)
        dr7filtered = self.filter_from_external(dr7_file,do_filter=False,return_id=True,keyword='gal_objID',exclude=False,mode='galaxy')
        dummy.cluster_galaxies = dr7filtered
        bossfiltered = self.filter_from_external(boss_file,do_filter=False,return_id=True,keyword='gal_objID',exclude=False,mode='galaxy')

        extra_data = {}
        if extras != None:
            for key in extras.keys():
                id = extras[key]['name']
                file = extras[key]['file']
                keyword = extras[key]['keyword']
                keyword2 = ''
                if extras.has_key('keyword2'):
                    keyword2 = extras[key]['keyword2']
                extra_filter = self.filter_from_external(file,do_filter=False,return_id=True,keyword=keyword,exclude=False,mode='galaxy',keyword2=keyword2)
                extra_data[id] = {'targets':extra_filter}


        gals = []
        avoid = ['total','BCG','new','available','detectable','addedBCG']
        if supplemental:
            avoid.append('grand_total')
            avoid.append('supplemental')
            
        for c in self.clusters:
            galaxies = c.galaxies
            if uber:
                galaxies = c.uber_cluster.galaxies
            c.boss_targets = {'total':0,'boss':[],'dr7':[],'2slaq':[],'wiggleZ':[],'new':[],'available':[],'detectable':[],'BCG':False,'addedBCG':False}
            if supplemental:
                c.boss_targets['supplemental'] = []
                c.boss_targets['grand_total'] = 0

            for key in extra_data.keys():
                c.boss_targets[key] = []
                
            register = []

            for g in galaxies:
                if g.objID in bossfiltered:
                    c.boss_targets['boss'].append(g)
                    c.boss_targets['total'] += 1
                elif g.objID in dr7filtered:
                    c.boss_targets['dr7'].append(g)
                    c.boss_targets['total'] += 1
                elif numpy.isnan(g.extras['2slaq']) == False:
                    c.boss_targets['2slaq'].append(g)
                    c.boss_targets['total'] += 1
                    
                elif numpy.isnan(g.extras['wiggleZ']) == False:
                    c.boss_targets['wiggleZ'].append(g)
                    c.boss_targets['total'] += 1

                elif (len(extra_data.keys()) > 0) and (True in [g.objID in extra_data[_id]['targets'] for _id in extra_data.keys()]):
                    keys = extra_data.keys()
                    test = [g.objID in extra_data[key]['targets'] for key in keys]
                    i = indexOf(test,True)
                    key = keys[i]

                    if g.objID in extra_data[key]['targets']:
                        if g.objID in register:
                            g.print_info()
                            g.parent.print_info()
                            return register
                        
                        c.boss_targets[key].append(g)
                        c.boss_targets['total'] += 1
                        register.append(g.objID)

                else:
                    c.boss_targets['available'].append(g)
                    if g.magnitudes[band] < maglim and g.magnitudes[band] > brightlim:
                        c.boss_targets['detectable'].append(g)
#                        if numpy.isnan(g.extras['2slaq']) and numpy.isnan(g.extras['wiggleZ']):
                    gals.append(g)


            for key in c.boss_targets.keys():
                if (key in avoid) == False:
                    if c.bcg in c.boss_targets[key]:
                        c.boss_targets['BCG'] = True


        if no_list:
            return register


        avoid = ['total','BCG','available','detectable','addedBCG']
        if supplemental:
            avoid.append('grand_total')
            avoid.append('supplemental')
        prioritised = []
        sup_list = []
        print 'Creating prioritisation list'
        from dataIO import Progress
        status = Progress(len(self.clusters))
        for c in self.clusters:
            priority = 0
            candidates = numpy.array(c.boss_targets['detectable'])
            if len(candidates) == 0:
                continue
#            candidates = numpy.array(c.boss_targets['available'])
            mag_test = numpy.array([g.magnitudes[band] < maglim and g.magnitudes[band] > brightlim for g in candidates])
            candidates = list(candidates[mag_test])
            candidates = list(candidates)
            if (c.bcg in candidates) and (c.bcg not in c.boss_targets['new']):
                c.boss_targets['new'].append(c.bcg)
                candidates.remove(c.bcg)
#                c.boss_targets['available'].remove(c.bcg)
                c.boss_targets['total'] += 1
                c.boss_targets['addedBCG'] = True
                c.bcg.priority = priority
                priority += 1
                c.bcg.sourcetype = 'BCG'

#            ncandidates = len(candidates)
            while (c.boss_targets['total'] < target) and len(candidates) > 0 and bcg_only == False:
 #               candidates = c.boss_targets['available']
                mags = {}
                for g in candidates:
                   mags[g.magnitudes[band]] = g
                choice = mags[min(mags.keys())]
                c.boss_targets['new'].append(choice)
#                c.boss_targets['available'].remove(choice)
                candidates.remove(choice)
                c.boss_targets['total'] += 1
                choice.priority = priority
                priority += 1
                choice.sourcetype = 'CGAL'
                if choice in self.uber_cluster_galaxies:
                    choice.sourcetype = 'ACGAL'

            prioritised = prioritised + list(c.boss_targets['new'])
            for key in c.boss_targets.keys():
                if (key in avoid) == False:
                    if c.bcg in c.boss_targets[key]:
                        c.boss_targets['BCG'] = True

            if supplemental:
                _keys = []
                for key in c.boss_targets.keys():
                    if (key in avoid) == False:
                        _keys.append(key)
                        
                band2 = band
                band2 = c.parent.selection['colours'][c.basis]['magnitude']
                supps = []
                uc = c.uber_cluster
                for g in uc.galaxies:
                    test = [g not in c.boss_targets[k] for k in _keys]
                    if (g not in c.galaxies) and countOf(test,True) == len(test):
                        if (g.magnitudes[band2] < c.bcg.magnitudes[band2]) and (g.magnitudes[band2] > brightlim):
                            supps.append(g)
                if len(supps) > 0:
                    smags = [g.magnitudes[band2] for g in supps]
                    g = supps[indexOf(smags,min(smags))]
                    c.boss_targets['supplemental'] = [g]
                    sup_list.append(g)
                c.boss_targets['grand_total'] = int(c.boss_targets['total'] + len(supps))
                
                
            status.update(0)
                    
                    
        
#        dummy = det(None,deepcopy(self.selection),dummy=True)
#        dr7filtered = self.filter_from_external(dr7_file,do_filter=True,keyword='gal_objID',exclude=True,mode='galaxy')
#        dummy.cluster_galaxies = dr7filtered
#        bossfiltered = dummy.filter_from_external(boss_file,do_filter=True,keyword='gal_objID',exclude=True,mode='galaxy')
#        bossfiltered = numpy.array(bossfiltered)
                    
        gals = numpy.array(gals)
        mags = numpy.array([g.magnitudes[band] for g in gals])
        filter = (mags < maglim) & (mags > brightlim)
        gals = gals[filter]
        print 'Updating galaxy target parameters'
        status = Progress(len(gals))
        for g in gals:
            try:
                pri = g.priority
            except:
                g.priority = 0
            if g == g.parent.bcg:
                g.sourcetype = 'BCG'
                g.cluster_bcg = True
            elif g in g.parent.associated_members:
                g.sourcetype = 'ACGAL'
                g.cluster_bcg = False
            elif g in sup_list:
                g.sourcetype = 'sBCG'
                g.cluster_bcg = False
                g.priority = 1
            else:
                g.sourcetype = 'CGAL'
                g.cluster_bcg = False
            status.update(0)
        for g in prioritised:
            if g in bcgs:
                g.cluster_bcg = True
            else:
                g.cluster_bcg = False

        for g in sup_list:
            g.sourcetype = 'sBCG'
            g.cluster_bcg = False
            g.priority = 1
                
            

        if write == True:        
            from tools import gallist_cat
            to_write = list(gals)+list(prioritised)
            if len(sup_list) > 0:
                to_write = to_write+list(sup_list)

            to_write = numpy.array(to_write)
            filter = numpy.array([numpy.isnan(g.camcol) for g in to_write])
            missing = to_write[filter]
            if len(missing) > 0:
                self.sdssimg_parameters(obj_override=list(missing),verbose=True)
            to_write = list(to_write)
            
            for g in to_write:
                g.id = g.objID
                g.magnitudes['fiber3mag_u'] = g.magnitudes['fiber_u']
                g.magnitudes['fiber3mag_g'] = g.magnitudes['fiber_g']
                g.magnitudes['fiber3mag_r'] = g.magnitudes['fiber_r']
                g.magnitudes['fiber3mag_i'] = g.magnitudes['fiber_i']
                g.magnitudes['fiber3mag_z'] = g.magnitudes['fiber_z']

                g.magnitudes['fiber2mag_u'] = 30.0 #g.magnitudes['fiber_u']
                g.magnitudes['fiber2mag_g'] = 30.0 #g.magnitudes['fiber_g']
                g.magnitudes['fiber2mag_r'] = 30.0 #g.magnitudes['fiber_r']
                g.magnitudes['fiber2mag_i'] = g.magnitudes['fiber_i2as']
                g.magnitudes['fiber2mag_z'] = 30.0 #g.magnitudes['fiber_z']

            
            mags = ['fiber2mag_u','fiber2mag_g','fiber2mag_r','fiber2mag_i','fiber2mag_z','fiber3mag_u','fiber3mag_g','fiber3mag_r','fiber3mag_i','fiber3mag_z']
            extras=['sourcetype','priority','run','rerun','camcol','field','id']
            exclude = []
            out2 = gallist_cat(gals,output='full_targets.fits',mags=mags,extras=extras,exclude=exclude,mag_suffix='')
            out2 = gallist_cat(prioritised,output='prioritised_targets.fits',mags=mags,extras=extras,exclude=exclude,mag_suffix='')
            if len(sup_list) > 0:
                out2 = gallist_cat(sup_list,output='supplementary_targets.fits',mags=mags,extras=extras,exclude=exclude,mag_suffix='')
                out2 = gallist_cat(list(prioritised)+list(sup_list),output='priority_supplemental.fits',mags=mags,extras=extras,exclude=exclude,mag_suffix='')
            for g in to_write:
                try:
                    del g.id
                    del g.magnitudes['fiber2mag_u']
                    del g.magnitudes['fiber2mag_g']
                    del g.magnitudes['fiber2mag_r']
                    del g.magnitudes['fiber2mag_i']
                    del g.magnitudes['fiber2mag_z']

                    del g.magnitudes['fiber3mag_u']
                    del g.magnitudes['fiber3mag_g']
                    del g.magnitudes['fiber3mag_r']
                    del g.magnitudes['fiber3mag_i']
                    del g.magnitudes['fiber3mag_z']
                except:
                    pass
        return gals,prioritised,sup_list

    def boss_targets(self,boss_file,keyword='gal_objID',fibre=1.0,uber=True,assoc=False,band='i',exclusion_list=None,max_req=2):
        from clusters import cluster,galaxies_within_r
        import numpy
        boss_matched = self.filter_from_external(boss_file,return_id=True,do_filter=False,keyword=keyword)
        self.compare('../boss/boss_only.cat',type='galaxies',colour='r',shape='diamond',name='boss_detected')
        boss_detected = self.external_data['boss_detected']['data']
        targets = []
        if exclusion_list != None:
            boss_matched = list(boss_matched) + list(exclusion_list)
            boss_matched = unique(boss_matched)
        for _c in self.clusters:
            cl_targ = []
            c = _c
            galaxies = c.galaxies
            if uber:
                try:
                    c = _c.uber_cluster
                except:
                    _c.generate_uber_cluster()
                    c = _c.uber_cluster
                galaxies = c.galaxies
            elif assoc:
                galaxies = c.galaxies+list(c.associated_members)

            filter = numpy.array([(g.objID in boss_matched) == False for g in galaxies])
            anti_filter = numpy.array([(g.objID in boss_matched) for g in galaxies])
            galaxies = numpy.array(galaxies)
            undet = galaxies[filter]
            bin = galaxies[anti_filter]

            boss_nearby = galaxies_within_r(c,boss_detected,r=c.max_radius()*4.0)


            exclusions = list(bin)+list(boss_nearby)
            cl_targets = []
#            cl_targets = list(bin)

            if len(undet) > 0:
                cl_tmp = cluster(undet,c.dna)
                gal_dict = {}

                end_state = False
                while len(cl_targets)+len(bin) < max_req:
#                while len(cl_targets) < max_req:
                    gal_dict = {}
                    for g in cl_tmp.galaxies:
                        if g not in exclusions:
                            gal_dict[g.magnitudes[band]] = g
                
                    mag = gal_dict.keys()
                    mag.sort()
                    cycle = iter(mag)
                    while 1:
                        try:
                            _g = gal_dict[cycle.next()]
                            test = [len(galaxies_within_r(avoid,[_g],r=radians(1.0/60.0))) == 0 for avoid in exclusions]
                            if (False in test) == False or len(test) == 0:
                                cl_targets.append(_g)
                                exclusions.append(_g)
                                break
                        except StopIteration:
                            end_state = True
                            break
                    if end_state:
                        break
            targets = targets + list(cl_targets)
            _c.re_posess()
        return targets
            
                    
    def wrap_restore():
        if self.prev_wrap_mode:
            self.wrap_toggle(mode=self.prev_wrap_mode)
        

    def wrap_toggle(self,mode=None):
        if mode in ['obs','observe']:
            mode = 'observed'
        elif mode in ['det','detect']:
            mode = 'detected'


        from math import radians
        try:
            self.prev_wrap_mode = self.wrap_mode
        except:
            self.prev_wrap_mode = self.clusters[0].wrap_mode
        try:
            pre_mode = self.clusters[0].wrap_mode
        except:
            pre_mode = None
        if self.clusters:
            for cl in self.clusters:
                cl.wrap_toggle(mode=mode)
        if self.external_data:
            for key in self.external_data.keys():
                for cl in self.external_data[key]['data']:
                    cl.wrap_toggle(mode=mode)

        self.wrap_mode = mode

        if self.selection:
            cx = [(cl.x < radians(360.0) and (cl.x >= 0.0)) for cl in self.clusters]
#            cx = [cl.x < radians(360.0) for cl in self.clusters]
            if False in cx:
                mode = 'detected'
                #need to -360
            else:
                mode = 'observed'
#            print mode

#            if mode == 'detected' and not True in [self.selection['limits']['xmin']>radians(360.0),self.selection['limits']['xmax']>radians(360.0)]:
            if mode == 'detected' and not True in [(self.selection['limits']['xmin']>radians(360.0) or self.selection['limits']['xmin'] < 0.0),(self.selection['limits']['xmax']>radians(360.0) or self.selection['limits']['xmin'] < 0.0)]:

                adjustment = radians(360.0)
                if self.selection['io']['dataset'] == 's82' and self.selection['limits']['xmin'] > 2.0:
#                if self.selection['io']['dataset'] == 's82' and self.selection['detection']['subset_index'] == -1 and self.selection['limits']['xmin'] > 2.0:
                    adjustment = -1.0*radians(360.0)


                print 'wrapping selection %f'%(degrees(adjustment))
                self.selection['limits']['xmin'] = self.selection['limits']['xmin'] + adjustment
                self.selection['limits']['xmax'] = self.selection['limits']['xmax'] + adjustment




#                print 'wrapping selection +360.0'
#                self.selection['limits']['xmin'] = self.selection['limits']['xmin'] + radians(360.0)
#                self.selection['limits']['xmax'] = self.selection['limits']['xmax'] + radians(360.0)

                if self.selection.has_key('grid'):
                    keys = ['xc','xmax','xmin']
                    convert_mode = lambda x: x
                    if self.selection['grid']['self']['units'] == 'radians':
                        from numpy import radians
                        convert_mode = radians

                    for key in keys:
                        self.selection['grid']['self'][key] = self.selection['grid']['self'][key] + convert_mode(360.0)
                    zones = ['always','choose','never']
                    try:
                        for zone in zones:
                            _zone = []
                            for sub_zone in self.selection['grid']['self']['zones'][zone]:
                                _sz = sub_zone
                                _sz[0] = _sz[0] + convert_mode(360.0)
                                _zone.append(_sz)
                            self.selection['grid']['self']['zones'][zone] = _zone
                                

                    except KeyError:
                        pass

                if self.selection.has_key('grid') and self.selection['grid']['self'].has_key('boundary'):
                    from numpy import radians
                    convert_mode = radians

                    master_boundary = []
                    for b in self.selection['grid']['self']['boundary']:
                        if type(b[0]) == type([]):
                            _b = []
                            for sub_bound in b:
                                b[0] = b[0] + convert_mode(360.0)
                                _b.append(b)
                            b = _b
                        else:
                            b[0] = b[0] + convert_mode(360.0)
                        master_boundary.append(b)
                       
                        
                    

            elif mode == 'observed' and True in [(self.selection['limits']['xmin']>radians(360.0) or self.selection['limits']['xmin'] < 0.0),(self.selection['limits']['xmax']>radians(360.0) or self.selection['limits']['xmax'] < 0.0)]:

                adjustment = -1.0*radians(360.0)
                if True in [self.selection['limits']['xmin'] < 0.0,self.selection['limits']['xmax'] < 0.0]:
                    adjustment = radians(360.0)

                print 'wrapping selection %f'%(degrees(adjustment))
                print self.selection['limits']['xmin']



                self.selection['limits']['xmin'] = self.selection['limits']['xmin'] + adjustment
                self.selection['limits']['xmax'] = self.selection['limits']['xmax'] + adjustment

                if self.selection.has_key('grid'):
                    keys = ['xc','xmax','xmin']
                    convert_mode = lambda x: x
                    if self.selection['grid']['self']['units'] == 'radians':
                        from numpy import radians
                        convert_mode = radians

                    for key in keys:
                        self.selection['grid']['self'][key] = self.selection['grid']['self'][key] - convert_mode(360.0)
                    zones = ['always','choose','never']
                    try:
                        for zone in zones:
                            _zone = []
                            for sub_zone in self.selection['grid']['self']['zones'][zone]:
                                _sz = sub_zone
                                _sz[0] = _sz[0] - convert_mode(360.0)
                                _zone.append(_sz)
                            self.selection['grid']['self']['zones'][zone] = _zone
                    except KeyError:
                        pass

                if self.selection.has_key('grid') and self.selection['grid']['self'].has_key('boundary'):
                    from numpy import radians
                    convert_mode = radians
                    master_boundary = []
                    for b in self.selection['grid']['self']['boundary']:
                        if type(b[0]) == type([]):
                            _b = []
                            for sub_bound in b:
                                b[0] = b[0] - convert_mode(360.0)
                                _b.append(b)
                            b = _b
                        else:
                            b[0] = b[0] - convert_mode(360.0)
                        master_boundary.append(b)
                        
                        
                        


                print
                print self.selection['limits']['xmin']
                print
                print 'Fin'


    def plot(self,perfect=True,ngals=5,galsize=4,label='',interact=False,cut=None,print_friendly=False,field=False,field_env=False,max_viz=False,presentation_mode=False,new_limits=None,ext_data=False,ingest_clusters=None,centroids=False,rings=False,colour='blue',boundary_colour='green',units='default',perfect_colour='!grey',observer_mode=True,plot_grid=True,source_background=False,fill=None,ra_buffer=10.0,lw=1.0,wrap_mode='None',bcg_colour='yellow',passback=False,auto_aspect=True,buttons=3):
        from detections import plot_clusters,plot_galaxies
        
        if passback == False:
            if auto_aspect:
                nmax = 19
                xmin = self.selection['limits']['xmin']
                xmax = self.selection['limits']['xmax']
                ymin = self.selection['limits']['ymin']
                ymax = self.selection['limits']['ymax']
                dx = xmax-xmin
                dy = ymax-ymin
                if (dx/dy > 0.8) and (dx/dy < 1.2):
                    fig = plt.figure()
                elif dx > dy:
                    asp = (nmax,nmax/(dx/dy)+2)
                    fig = plt.figure(figsize=asp)
                else:
                    asp = (nmax/(dx/dy)+2,nmax)
                    fig = plt.figure(figsize=asp)
            else:
                fig = plt.figure()

        if wrap_mode != 'None':
            self.wrap_toggle(mode=wrap_mode)


#        print [cl.x for cl in self.external_data.values()[0]['data']]

        if presentation_mode:
            size = 0.0
            bg = 'white'
            boundary_colour='black'
            units = 'degrees'
            if passback == False:sp = plt.subplot(111, axisbg=bg,aspect='equal')

            plot_bound = False
            try:
                bound = self.selection['boundary']
                plot_bound = True
            except:
                pass

            if (self.selection['detection'].has_key('boundary') or plot_bound) and (source_background):
                if passback == False:sp = plt.subplot(111, axisbg='white',aspect='equal')
                self.plot_boundary(ec='green',ls='solid',lw=2.0)

            elif self.selection['detection'].has_key('boundary') or plot_bound:
                if passback == False:sp = plt.subplot(111, axisbg='white',aspect='equal')
                self.plot_boundary(ec='green',ls='solid',lw=2.0,fc='None')
            else:
                if passback == False:sp = plt.subplot(111, axisbg=bg,aspect='equal')

            if self.selection.has_key('grid') and plot_grid:
                if not self.selection['grid'].has_key('no plot'):
                    try:
                        self.plot_grid(colour='green')
                    except KeyError:
                        pass

            limits = self.selection
            if new_limits != None:
                limits = new_limits
            
            if ext_data or self.external_data != None:
                self.extdata_plot(passback=True,colour='r',points=False,lw=lw)

            radius = None
            if perfect and 'mock' in self.selection['io']['dataset']:
                radius = max([abs(self.selection['limits']['xmin']),abs(self.selection['limits']['xmax']),abs(self.selection['limits']['ymin']),abs(self.selection['limits']['ymax'])])
                try:
                    self.selection['overrides']['no_fd'] = True
                    self.selection['overrides']['no_vir'] = True
                    
                except:
                    self.selection['overrides'] = {}
                    self.selection['overrides']['no_fd'] = True
                    self.selection['overrides']['no_vir'] = True
    
            if len(self.field_galaxies) > 0:
                plot_galaxies(self.field_galaxies,limits=limits,passback=True,colour='gray',size=size,verbose=True,fill='white',lw=lw)
            bcgs = [c.bcg for c in self.clusters]

            #plot in order of area:
            from science import science
            sci = science()
            filter = sci.filter(self,'area',0,9999,rank=True,reverse=True)

            plot_clusters(filter.clusters,colour='blue',limits=limits,passback=True,size=size,field_gals=None,max_viz=False,fill='blue',ext=False)
            if bcg_colour != None:
                plot_galaxies(bcgs,limits=limits,passback=True,colour=bcg_colour,size=size,verbose=False,fill=bcg_colour)

            if perfect and 'mock' in self.selection['io']['dataset']:
                plot_clusters(self.perfect_clusters,rings=True,centroids=True,passback=True,colour=perfect_colour)
                del self.selection['overrides']['no_fd']
                del self.selection['overrides']['no_vir']

            if radius != None:
                if observer_mode:
                    plt.xlim(radius,-1.0*radius)
                else:
                    plt.xlim(-1.0*radius,radius)
            
                plt.ylim(-1.0*radius,radius)
                cir = plt.Circle((0.0,0.0), radius=radius,fc='None',alpha=1.0,lw=1.0,ec='k')
#                cir = plt.Circle((0.0,0.0), radius=radius,fc='None',alpha=1.0,lw=lw,ec='k')
                plt.gca().add_patch(cir)


            elif observer_mode:
                xmin,xmax = plt.xlim()
                plt.xlim(xmax,xmin)


            xlim = plt.xlim()
            ylim = plt.ylim()
            

            if units != 'default':
                from numpy import radians
                if units == 'degrees':
                    xmin,xmax = plt.xlim()
                    ymin,ymax = plt.ylim()

                    xmin,xmax=degrees(xmin),degrees(xmax)
                    dr = float(xmax-xmin)/7.0

                    x_radians = plt.xticks()[0]
                    xvals = numpy.array([degrees(x_r) for x_r in x_radians])

                    y_radians = plt.yticks()[0]
                    yvals = [degrees(y_r) for y_r in y_radians]

                    if xvals[-1] < 0:
                        from math import pi
                        xvals = xvals + (360.0)
                    newx = ['%2.1f'%(_x) for _x in xvals]

                    newy = ['%2.1f'%(_y) for _y in yvals]
                    plt.xticks(x_radians,newx,fontsize=16)
                    plt.yticks(y_radians,newy,fontsize=16)
                   

            else:
                units = 'radians'
                
        

            plt.xlim(xlim)
            plt.ylim(ylim)

            if ra_buffer:
                from math import radians
                plt.xlim(self.selection['limits']['xmin']-radians(ra_buffer/60.0),self.selection['limits']['xmax']+radians(ra_buffer/60.0))
                plt.ylim(self.selection['limits']['ymin']-radians(ra_buffer/60.0),self.selection['limits']['ymax']+radians(ra_buffer/60.0))


                
            if observer_mode:
                xmin,xmax = plt.xlim()
                plt.xlim(max([xmax,xmin]),min([xmax,xmin]))

            plt.subplots_adjust(top=0.99,bottom=0.05,left=0.12,right=1.00)
            plt.xlabel('RA (%s)'%(units),fontsize=16)
            plt.ylabel('DEC (%s)'%(units),fontsize=16)
            
            return    

        if self.selection['io']['dataset'] == 'reconstructed':
            rings = True
            centroids = True
        bg = 'black'

        if print_friendly:
            bg = 'white'
            boundary_colour='green'

        if source_background:
            bg = 'white'
            boundary_colour='green'


        sp = plt.subplot(111, axisbg=bg,aspect='equal')

        if source_background:
            try:
                data = self.cuts['x'][0]
                data = self.cuts
            except:
                pass
                from copy import deepcopy
                tmp_sel = deepcopy(self.selection['colours'])
                for key in self.selection['colours'].keys():
                    self.selection['colours'][key]['width'] = 100.0
                self.get_input()
                data = self.cuts
                self.selection['colours'] = tmp_sel
            plt.hexbin(data['spatial']['x'],data['spatial']['y'],gridsize=200,cmap=plt.cm.binary,antialiased=True)

######            plt.hexbin(data['spatial']['x'],data['spatial']['y'],gridsize=200,bins='log',cmap=plt.cm.binary,antialiased=True)
            fill = True
            galsize = 0.0

        plot_bound = False
        try:
            bound = self.selection['boundary']
            plot_bound = True
        except:
            pass

        if self.selection.has_key('boundary') or plot_bound:
            self.plot_boundary(ec=boundary_colour,fc='None',ls='solid',lw=2.0)
            
        if self.selection.has_key('grid') and plot_grid:
            if 1:
#            try:
                self.plot_grid(colour='white')
            else:
#            except:
                pass
        pg_x,pg_y = [],[]
        if perfect and 'mock' in self.selection['io']['dataset'] and len(self.perfect_clusters) > 0:
            plot_clusters(self.perfect_clusters,rings=True,centroids=True,passback=True,colour=perfect_colour)
            
        if len(self.field_galaxies) > 0 and field:
            print '\tPlotting field_galaxies - this may take some time:'
            plot_galaxies(self.field_galaxies,limits=self.selection,passback=True,colour='gray',size=galsize,verbose=True)

        field_gals = None
        if self.field_galaxies and field_env:
            field_gals = self.field_galaxies

        if (len(self.field_galaxies) > 0 and field) or field_env == False:
            #will have already plotted these 
            field_gals = None
            
            
        plot_clusters(self.clusters,colour=colour,limits=self.selection,passback=True,size=galsize,field_gals=field_gals,max_viz=max_viz,centroids=centroids,rings=rings,ext=False,fill=fill,bcg_colour=bcg_colour)
        
        try:
            if self.external_data != None:
                self.extdata_plot(passback=True,colour='r',lw=lw)
        except:
            pass
        if label == '':
            if self.selection['overrides'].has_key('rho_crit'):
                rc = self.selection['overrides']['rho_crit']
            else:
                try:
                    rc = self.selection['detection']['rho_crit']
                except:
                    rc = -99.0
                    
            if self.selection['overrides'].has_key('probthresh'):
                pt = self.selection['overrides']['probthresh']
            else:
                try:
                    pt = self.selection['detection']['probthresh']
                except:
                    pt = -99.0
            
            try:
                label = '%s %3.3f sq. deg, P=%3.3f,rho_c=%3.3f'%(self.selection['defaults']['colour'],self.selection['detection']['detection_area'],pt,rc)
            except:
                label = ''

        plt.xlabel('RA')
        plt.ylabel('DEC')
        plt.title(label)
        plt.xlim(self.selection['limits']['xmin'],self.selection['limits']['xmax'])
        plt.ylim(self.selection['limits']['ymin'],self.selection['limits']['ymax'])
        print self.selection['limits']['xmin'],self.selection['limits']['xmax']
        if ra_buffer:
            from math import radians
            plt.xlim(self.selection['limits']['xmin']-radians(ra_buffer/60.0),self.selection['limits']['xmax']+radians(ra_buffer/60.0))
            plt.ylim(self.selection['limits']['ymin']-radians(ra_buffer/60.0),self.selection['limits']['ymax']+radians(ra_buffer/60.0))
        
        if observer_mode:
            xmin,xmax = plt.xlim()
            plt.xlim(xmax,xmin)

        plt.subplots_adjust(left=0.12,right=0.97,top=0.91,bottom=0.08)


        if units != 'default':
            xmin,xmax = plt.xlim()
            ymin,ymax = plt.ylim()

            from numpy import radians
            if units == 'degrees':
#                xmin,xmax=degrees(xmin),degrees(xmax)
#                ymin,ymax=degrees(ymin),degrees(ymax)

#                xmin,xmax=degrees(xmin),degrees(xmax)
#                ymin,ymax=degrees(ymin),degrees(ymax)
                
                dr = float(xmax-xmin)/7.0

                x_radians = plt.xticks()[0]
                xvals = numpy.array([float(degrees(x_r)) for x_r in x_radians])

                y_radians = plt.yticks()[0]
                yvals = [float(degrees(y_r)) for y_r in y_radians]

                if xvals[-1] < 0:
                    from math import pi
                    xvals = xvals + (360.0)
                newx = ['%2.1f'%(_x) for _x in xvals]

                newy = ['%2.1f'%(_y) for _y in yvals]
                plt.xticks(x_radians,newx)
                plt.yticks(y_radians,newy)

#                plt.xticks(x_radians,xvals)
 #               plt.yticks(y_radians,yvals)



            else:
                units = 'radians'
            print xmin,xmax
            print ymin,ymax
            plt.xlim(xmin,xmax)
            plt.ylim(ymin,ymax)


        if interact:
            action_button = 2
            if buttons == 2:
                action_button = 1
            plt.ioff()
#            plt.ion()
            plt.subplots_adjust(left=0.07,right=0.98,top=0.93,bottom=0.07)
#            sp.axes.set_position([-0.05,0.1,0.8,0.8])

            ii = interaction(fig,self,action_button=action_button)
#            plt.show()
            plt.ioff()
            plt.close('all')


    def plot_grid(self,grid='default',colour='w',lw=0.5,alpha=0.2):
#        print 'plottling grid'
        if grid == 'default':
            try:
                _this_grid = self.selection['grid']['self']
            except:
                print 'Could not plot grid - please specify one'
        else:
            _this_grid = grid

        if _this_grid.has_key('stripe'):
            convert_mode = lambda x: x
            if _this_grid['units'] == 'degrees':
                from numpy import radians
                convert_mode = radians
            colours = {'always':'g','choose':'y','never':'r'}
            for zone_type in _this_grid['zones']:
                if not zone_type == 'chooser':
                    for zone in _this_grid['zones'][zone_type]:
                        x0 = convert_mode(zone[0])
                        y0 = convert_mode(zone[1])
                        dx = convert_mode(zone[2])
                        dy = convert_mode(zone[3])
                        rec = plt.Rectangle((x0,y0),dx,dy,fc='None',lw=lw,ec=colours[zone_type],ls='dashed',alpha=alpha)
                        plt.gca().add_patch(rec)

                        

        return

        dr = _this_grid['dr']
        ddr = _this_grid['ddr']
        ddr = ddr/3.0
        xmin,ymin = _this_grid['xc']-_this_grid['dr'],_this_grid['yc']-_this_grid['dr']
        

        if _this_grid['units'] == 'degrees':
            import numpy
            xmin = numpy.radians(xmin)
            ymin = numpy.radians(ymin)
            dr = numpy.radians(dr)
            ddr = numpy.radians(ddr)
        xmin = xmin - ddr
        ymin = ymin - ddr
        rec = plt.Rectangle((xmin,ymin),(dr+ddr)*2,(dr+ddr)*2,fc='None',lw=lw,ec='g',ls='dashed',alpha=alpha)
        plt.gca().add_patch(rec)
        xmin = xmin - ddr
        ymin = ymin - ddr
        rec = plt.Rectangle((xmin,ymin),(dr+(2.0*ddr))*2,(dr+(2.0*ddr))*2,fc='None',lw=lw,ec='y',ls='dashed',alpha=alpha)
        plt.gca().add_patch(rec)
    def plot_boundary(self,alpha=0.25,lw=0.0,ec='grey',fc='gray',ls=None):

        if self.selection.has_key('boundary'):
            _col = '#D4D4D4'
            if type(self.selection['boundary']) == type({}):
                bb = self.selection['boundary']
                for key in bb.keys():
                    for patch in bb[key]:
                        p = plt.Rectangle(patch[0],patch[1],patch[2],fc=_col,alpha=1.0,lw=lw,ec='none',ls=ls)
#                        plt.
                        gca().add_patch(p)

            elif type(self.selection['boundary'][0][0]) == type([]) or type(self.selection['boundary'][0][0]) == numpy.ndarray:
                print 'This channel'

                boundary = self.selection['boundary']
                if type(boundary[0][0]) == type([]) or type(boundary[0][0]) == numpy.ndarray:
                    print 'Nest fix - change the source grid data to clear this!!!'
#                    boundary = boundary[0]
                for p in boundary:
                    patch = plt.Polygon(p,fc=fc,alpha=alpha,lw=lw,ec=ec,ls=ls)
#                    patch = plt.Polygon(p,alpha=alpha,lw=lw,ec=ec,fc=fc,ls=ls)
                    plt.gca().add_patch(patch)
            else:
                patch = plt.Polygon(self.selection['boundary'],fc=fc,alpha=alpha,lw=lw,ec=ec,ls=ls)
#                patch = plt.Polygon(self.selection['boundary'],alpha=alpha,lw=lw,ec=ec,fc=fc,ls=ls)
                plt.gca().add_patch(patch)

                
    def get_subids(self,grid,verbose=False):
        """\
        for all clusters in detection class, assign sub_id based on limits from gridfile

        """
        #set all cluster ids to None
        from dataIO import Progress
        import numpy
        for cl in self.clusters:
            cl.sub_id = None
        #which way? per cluster, or per region?
        #if len(self.clusters) > len(grid.keys()):
        status = Progress(len(grid.keys()))
        for key in grid.keys():
            sub_id = grid[key]['sub_id']
            dummy = self.clone()
            dummy.snap_to_grid(grid,id=key)
            dummy.clusters = self.clusters
            dummy.clusters = dummy.enforce_limits(add=False)
            if len(dummy.clusters) > 0:
                test = numpy.array([cl.sub_id == None for cl in dummy.clusters])
                clusters = numpy.array(dummy.clusters)
                if False in test:
                    pre_exist = clusters[test == False]
                    new = clusters[test]
                    for cl in pre_exist:
                        if type(cl.sub_id) == type([]):
                            cl.sub_id.append(sub_id)
                        else:
                            cl.sub_id = [cl.sub_id]
                            cl.sub_id.append(sub_id)
                    for cl in new:
                        cl.sub_id = sub_id
                else:
                    for cl in clusters:
                        cl.sub_id = sub_id
            status.update(0)
                     
        test = numpy.array([cl.sub_id == None for cl in self.clusters])
        if True in test:
            clusters = numpy.array(self.clusters)
            failed = clusters[test]
            print 'There are %d clusters that we did not get sub_ids for...'%(len(failed))
        else:
            print 'All done'
            


    def enforce_limits(self,verbose=False,add=True):
        """\
        Boots out any clusters assigned to this detection class that live outside of the defined survey limits

        """

        cx = numpy.array([c.x for c in self.clusters])
        cy = numpy.array([c.y for c in self.clusters])
        clusters = numpy.array(self.clusters)
        filterx = (cx <= self.selection['limits']['xmax']) & (cx >= self.selection['limits']['xmin'])
        filtery = (cy <= self.selection['limits']['ymax']) & (cy >= self.selection['limits']['ymin'])
        filter = filterx & filtery
        if verbose:
            print 'Reduced clusters from %d to %d in this detection class'%(len(clusters),countOf(filter,True))
        clusters = list(clusters[filter])
        if add:
            self.addclusters(clusters)
        else:
            return clusters

        

    def snap_to_clusters(self,clusters=None,r='auto'):
        """\
        applies the characteristics/limits of either supplied clusters or hosted clusters to this detection class.
        r will either be cluster.max_radius+10percent, or a defined radius (in radians)

        """
        xmaxes = []
        ymaxes = []
        xmins = []
        ymins = []    
        for cl in self.clusters:
            cx = cl.x
            cy = cl.y
            max_radius = 1.1*cl.max_radius()
            if r != 'auto':
                max_radius = r
            xmaxes.append(cx+max_radius)
            ymaxes.append(cy+max_radius)
            xmins.append(cx-max_radius)
            ymins.append(cy-max_radius)


        xmin = min(xmins)
        ymin = min(ymins)
        xmax = max(xmaxes)
        ymax = max(ymaxes)

        self.selection['limits']['xmin'] = xmin
        self.selection['limits']['xmax'] = xmax
        self.selection['limits']['ymin'] = ymin
        self.selection['limits']['ymax'] = ymax
            


    def snap_to_cluster(self,cluster,r='auto'):
        """\
        applies the characteristics/limits of the supplied cluster to this detection class.
        r will either be cluster.max_radius+10percent, or a defined radius (in radians)

        """
        cl = cluster
        cx = cl.x
        cy = cl.y
        max_radius = 1.1*cl.max_radius()
        if r != 'auto':
            max_radius = r
        self.selection['limits']['xmin'] = cx-max_radius
        self.selection['limits']['xmax'] = cx+max_radius
        self.selection['limits']['ymin'] = cy-max_radius
        self.selection['limits']['ymax'] = cy+max_radius
        return
        
    def snap_to_selection(self,selection,clip=False):
        """\
        applies the limits of the supplied selection dictionary (eg from another detection) to this detection
        clip         : include galaxies only within these limits - remove others
        """
        self.selection['limits']['xmin'] = selection['limits']['xmin']
        self.selection['limits']['xmax'] = selection['limits']['xmax']
        self.selection['limits']['ymin'] = selection['limits']['ymin']
        self.selection['limits']['ymax'] = selection['limits']['ymax']
        if clip:
            inside = []
            cg = numpy.array(self.cluster_galaxies)
            cg_x = numpy.array([g.x for g in cg])
            cg_y = numpy.array([g.y for g in cg])
            filterx = (cg_x >= self.selection['limits']['xmin']) & (cg_x <= self.selection['limits']['xmax'])
            filtery = (cg_y >= self.selection['limits']['ymin']) & (cg_y <= self.selection['limits']['ymax'])
            filter = filterx & filtery
            cg = cg[filter]
            self.cluster_galaxies = list(cg)

    def expand_limits(self,dx,dy=None,units='arcmins'):
        """\
        applies the characteristics/limits of the supplied grid to this detection class

        """
        from numpy import radians
        if (dx == None) and (dy == None):
            print 'No changes requested' 
            return

        conv = {}
        conv['arcmins'] = 60.0
        conv['arcsecs'] = 3600.0
        conv['degrees'] = 1.0


        dx_r = radians(dx/conv[units])
        dy_r = radians(dx/conv[units])
        if dy:
            dy_r = radians(dy/conv[units])
        self.selection['limits']['xmin'] -= dx_r
        self.selection['limits']['xmax'] += dx_r
        self.selection['limits']['ymin'] -= dy_r
        self.selection['limits']['ymax'] += dx_r
        return
        
    def plot_limits(self):
        """\
        sets plot x/y limits to be the survey limits.

        """
        import matplotlib.pyplot as plt
        plt.xlim(self.selection['limits']['xmax'],self.selection['limits']['xmin'])
        plt.ylim(self.selection['limits']['ymin'],self.selection['limits']['ymax'])

        

    def snap_to_grid(self,grid,id=None):
        """\
        applies the characteristics/limits of the supplied grid to this detection class

        """
        source_grid = grid
        if id:
            try:
                grid = source_grid[id]
            except:
                if type(source_grid.values()[0]) == type({}):
                    print 'Searching through stripes...'
                    bail = False
                    for key in source_grid.keys():
                        if id in source_grid[key].keys():
                            src_grid = source_grid[key]
                            _grid = source_grid[key][id]
                    try:
                        source_grid = src_grid
                        grid = _grid
                    except:
                        raise SystemExit('#snap_to_grid(): Cannot find id %d'%(id))
                else:
                    raise SystemExit('#snap_to_grid(): Cannot find id %d'%(id))
        d = self
        d.selection['grid'] = {}
        try:
            d.selection['detection']['subset_index'] = grid['sub_id']
        except:
            pass
        if source_grid != grid:
            d.selection['grid'] = deepcopy(source_grid)
        d.selection['grid']['self'] = grid
        try:
            d.selection['detection']['stripe'] = grid['stripe']
        except:
            pass
        
        d.selection['detection']['subset_index'] = grid['sub_id']
        d.set_limits()
        if not d.clusters:
            d.clusters = []
    
        return

    def set_limits(self,grid='default',area=True):
        from math import pi
        if grid == 'default':
            try:
                _this_grid = self.selection['grid']['self']
            except:
                print 'Please specifiy a grid'
        else:
            _this_grid = grid
        xmin,xmax = _this_grid['xmin'],_this_grid['xmax']
        ymin,ymax = _this_grid['ymin'],_this_grid['ymax']
        if _this_grid['units'] == 'degrees':
            xmin = numpy.radians(xmin)
            xmax = numpy.radians(xmax)
            ymin = numpy.radians(ymin)
            ymax = numpy.radians(ymax)
        self.selection['limits']['xmin'],self.selection['limits']['xmax'] = xmin,xmax
        self.selection['limits']['ymin'],self.selection['limits']['ymax'] = ymin,ymax
        if area:
            self.selection['detection']['detection_area_radians'] = (xmax-xmin)*(ymax-ymin)
            self.selection['detection']['detection_area'] = (xmax-xmin)*(ymax-ymin)*((360.0*360.0)/(4*pi*pi))

            if _this_grid.has_key('area_radians'):
                self.selection['overrides']['area'] = _this_grid['area_radians']
                self.selection['detection']['detection_area_radians'] = _this_grid['area_radians']
                if _this_grid.has_key('area_degrees'):
                    self.selection['detection']['detection_area'] = _this_grid['area_degrees']
                else:
                    from math import pi
                    self.selection['detection']['detection_area'] = _this_grid['area_radians']


        if _this_grid.has_key('boundary'):
            self.selection['boundary'] = _this_grid['boundary']

        try:
            from numpy import degrees
            from clusters import clusters_within_r
            ed = self.external_data
            dx,dy = self.selection['limits']['xmax']-self.selection['limits']['xmin'],self.selection['limits']['ymax']-self.selection['limits']['ymin']
            xc,yc = self.selection['limits']['xmin']+(dx/2.0),self.selection['limits']['ymin']+(dy/2.0)
            if type(ed) == type({}):
                for key in ed.keys():
                    src = numpy.array(ed[key]['data'])
                    if self.selection.has_key('grid'):
                        ra = numpy.array(ed[key]['x'])
                        dec = numpy.array(ed[key]['y'])
                        filter = self.boundary_check(ra,dec)
                        nearby = src[filter]

                    else:
                        nearby = clusters_within_r([xc,yc],src,r=dx,shape='square')

                    ed[key]['data'] = numpy.array(nearby)
                    if ed[key].has_key('r'):
                        ed[key]['r'] = numpy.array([degrees(cl.extent)*60.0 for cl in ed[key]['data']])
                    ed[key]['x'] = numpy.array([cl.x for cl in ed[key]['data']])
                    ed[key]['y'] = numpy.array([cl.y for cl in ed[key]['data']])
                    
        except:
            print 'Warning: could not set detection selection limits from grid data'

    def merge(self,sub_ids='default',keep_ownership=True):
        from dataIO import load
        from tools import reconstruct
        if sub_ids == 'default':
            try:
                sub_ids = self.selection['grid']['self']['neighbourN']
            except:
                print 'Cannot determine neighbours - aborting'
                return
        grid = load('mygrid.dat')
        print 'Need to use a better grid!!!'

        det = reconstruct('../dr7/dr7full_test.fit',dataset='dr7full',grid=grid)
        det.compare('../maxbcg.cat',colour='red',shape='circle',points=False,type='cluster',name='maxBCG')
        det.compare('../abell.dat',colour='yellow',shape='square',type='cluster',name='Abell',points=False)
        neighbours = [det.reconstruct(id) for id in sub_ids]
        clusters = self.clusters
        boundary = []
        if self.selection.has_key('boundary'):
            if type(self.selection['boundary'][0][0]) == type([]):
                boundary = boundary + self.selection['boundary']
            else:
                boundary.append(self.selection['boundary'])

        for d in neighbours:
            clusters = clusters + d.clusters
            if d.selection.has_key('boundary'):
                if type(d.selection['boundary'][0][0]) == type([]):
                    boundary = boundary + d.selection['boundary']
                else:
                    boundary.append(d.selection['boundary'])
            
        output = self.clone()
        output.mirror = det.mirror
        
        output.selection['grid']['self']['xmin'] = min([d.selection['grid']['self']['xmin'] for d in neighbours])
        output.selection['grid']['self']['xmax'] = max([d.selection['grid']['self']['xmax'] for d in neighbours])
        output.selection['grid']['self']['ymin'] = min([d.selection['grid']['self']['ymin'] for d in neighbours])
        output.selection['grid']['self']['ymax'] = max([d.selection['grid']['self']['ymax'] for d in neighbours])
        try:
            output.selection['grid']['self']['area_radians'] = output.selection['grid']['self']['area_radians'] + sum([d.selection['grid']['self']['area_radians'] for d in neighbours])
        except:
            pass
        try:
            output.selection['grid']['self']['area_degrees'] = output.selection['grid']['self']['area_degrees'] + sum([d.selection['grid']['self']['area_degrees'] for d in neighbours])
        except:
            pass
        output.selection['detection']['subset_indices'] = [self.selection['detection']['subset_index']] + sub_ids
        output.external_data = det.external_data
        output.set_limits()
        output.addclusters(clusters)
        if len(boundary) > 0:
            output.selection['boundary'] = boundary

        if keep_ownership:
            for d in neighbours:
                d.re_posess()
        
            
        return output

        #grids

        #limits

        #boundaries
            
    
    

    def reconstruct(self,sub_id,vertices=True,verbose=False,dataset='dr7full',clusterID=None):
        sid = sub_id
        if (sid == None) and (clusterID != None):
            from tools import dr7
            det = dr7()
            clusters = {}
            for cl in det.clusters:
                clusters[cl.dna] = cl
            if clusterID in clusters.keys():
                cl = clusters[clusterID]
                sid = cl.sub_id
                print 'Cluster is in sub_id %d'%(sid)
            else:
                print 'ClusterID %d not found in cluster list and no sub_id defined: aborting'%(clusterID)
                return None
                

        try:
            d = self.mirror[sid]
            if clusterID == None:
                for cl in d.clusters:
                    clusters[cl.dna] = cl
                if clusterID in clusters.keys():
                    d = clusters[clusterID]

        except:
            from dataIO import load
            path = '/export/disk/dmurphy/cluster_cats/%s/'%(dataset)
            if vertices:
                path = path+'full/'
            else:
                path = path+'no_vertex/'
            import os
            files = os.listdir(path)

            options = []
            for file in files:
                if '.dat' in file and '_region%d.'%(sid) in file:
                    options.append(path+file)
            if len(options) > 1:
                print 'More than one option to chose from'
                return None
            elif len(options) == 0:
                print 'No files to choose from'
                return None
            file = options[0]
            d = load(file)
            try:
                self.mirror[sid] = d
            except:
                print 'Creating .mirror link'
                self.mirror = {}
                self.mirror[sid] = d

            if self.selection.has_key('grid'):
                grid = self.selection['grid']
                keys = grid.keys()
                values = grid.values()
                sids = [grid[k]['sub_id'] for k in keys]
                if sub_id in sids:
                    j = indexOf(sids,sub_id)
                    d.selection['grid'] = {}
                    d.selection['grid']['self'] = values[j]
                    if d.selection['grid']['self'].has_key('boundary'):
                        d.selection['boundary'] = d.selection['grid']['self']['boundary']
                    
                    for k in d.selection['grid']['self']['neighbourN']: d.selection['grid'][k] = grid[keys[indexOf(sids,k)]]
                    d.set_limits(area=True)

        if clusterID != None:
            clusters = {}
            for cl in d.clusters:
                clusters[cl.dna] = cl
            if clusterID in clusters.keys():
                d = clusters[clusterID]
        return d


    def zdist(self,plot=True,bins='auto',label=None,colour='k',ls='-'):
        import scipy
        cz = [cl.z for cl in self.clusters]
        if bins == 'auto':
            bins = 20

        try:
            Nz,z = scipy.histogram(cz,new=False,bins=bins)
        except:
            Nz,z = scipy.histogram(cz,bins=bins)
            Nz = list(Nz)+[0]
        plt.plot(z,Nz,color=colour,linestyle=ls,label=label)
        plt.xlabel('z')
        plt.ylabel('N(z)')
        return [Nz,z]



    def zdist_old(self,plot=False,type=['clusters'],dz=0.05,zmax=1.55):
        source_gals = self.cuts['halo_zcos']
        cluster_gals = [galaxy.zcos for galaxy in self.cluster_galaxies]
        clusters = [cluster.zcos for cluster in self.clusters]
        perfect_clusters = [cluster.zcos for cluster in self.perfect_clusters]
        z = arange(0.0,zmax+(dz/2.0),dz)
        N_sg = zeros(len(z))
        N_cg = zeros(len(z))
        N_c = zeros(len(z))
        N_pc = zeros(len(z))
        for i in xrange(len(z)):
            lower = [galz > z[i] for galz in source_gals]
            upper = [galz <= z[i]+dz for galz in source_gals]
            test = [lower[j]*upper[j] == 1 for j in xrange(len(lower))]
            N_sg[i] = countOf(test,True)

            lower = [galz > z[i] for galz in cluster_gals]
            upper = [galz <= z[i]+dz for galz in cluster_gals]
            test = [lower[j]*upper[j] == 1 for j in xrange(len(lower))]
            N_cg[i] = countOf(test,True)

            lower = [clusterz > z[i] for clusterz in clusters]
            upper = [clusterz <= z[i]+dz for clusterz in clusters]
            test = [lower[j]*upper[j] == 1 for j in xrange(len(lower))]
            N_c[i] = countOf(test,True)

            lower = [clusterz > z[i] for clusterz in perfect_clusters]
            upper = [clusterz <= z[i]+dz for clusterz in perfect_clusters]
            test = [lower[j]*upper[j] == 1 for j in xrange(len(lower))]
            N_pc[i] = countOf(test,True)

        counts = {'clusters':N_c,'clustergalxies':N_cg,'galaxies':N_sg,'perfectclusters':N_pc}
        if len(unique(counts.keys()+type)) != len(counts.keys()):
            for t in type:
                if not t in counts.keys():
                    print '\t%s does not exist in this data'%(t)
                    return counts
                    
        if plot:
            plt.xlabel('redshift')
            plt.ylabel('N')
            if len(type) > 0: col = ['k','r','g','b']
            else : col = ['k']
            colours = iter(col)
            for t in type:
                plt.plot(z,counts[t],'-%s'%(colours.next()),label=t)
            plt.xlim(0.0,1.5)
            plt.legend()
            plt.show()
        return counts
        
    def perfect_outliers(self,ggroup,outlier_metric='r80',resample='bootstrap',nsigma=3.0,verbose=False,plot=False,no_clobber=False):
        from tools import bootstrap as resampler
        from clusters import cluster
        
        if self.selection['mock_data'].has_key('outlier_metric'):
            outlier_metric = self.selection['mock_data']['outlier_metric']

        if self.selection['mock_data'].has_key('outlier_significance'):
            nsigma = self.selection['mock_data']['outlier_significance']

        if self.selection['mock_data'].has_key('no_clobber'):
            no_clobber = self.selection['mock_data']['no_clobber']

        if self.selection['mock_data'].has_key('outlier_verbosity'):
            verbose = self.selection['mock_data']['outlier_verbosity']
        

        if resample == 'jacknife':
            from tools import jacknife as resampler
        include = False
        if 1:
            if 1:
                if outlier_metric == 'r80':
                    proto_cluster = cluster(ggroup,-1)
                    r80 = proto_cluster.r_x()
                    test = numpy.array([proto_cluster.distance_from(xc=g.x,yc=g.y) < 2.0*r80 for g in proto_cluster.galaxies])
                    if countOf(test,True) >=self.selection['detection']['smallest_cluster']:
                        proto_galaxies = numpy.array(proto_cluster.galaxies)
                        ggroup = proto_galaxies[test]
                        include = True

                if outlier_metric == 'virial':
                    from tools import ang_diam_dist
                    proto_cluster = cluster(ggroup,-1)
                    rvir = proto_cluster.r_vir
                    test = numpy.array([abs(ang_diam_dist(proto_cluster,g)) < rvir for g in proto_cluster.galaxies])
                    if countOf(test,True) >=self.selection['detection']['smallest_cluster']:
                        proto_galaxies = numpy.array(proto_cluster.galaxies)
                        ggroup = proto_galaxies[test]
                        include = True


                elif outlier_metric == 'radius':
                    proto_cluster = cluster(ggroup,-1)
                    r80 = proto_cluster.r_x()
#                    distances = numpy.array([proto_cluster.distance_from(xc=g.x,yc=g.y) < 2.0*r80 for g in proto_cluster.galaxies])
                    distances = numpy.array([proto_cluster.distance_from(xc=g.x,yc=g.y) for g in proto_cluster.galaxies])
                    dist_median,dist_sigma = resampler(distances)
                    test = numpy.array([d < dist_median+(nsigma*dist_sigma) for d in distances])
                    if countOf(test,True) >=self.selection['detection']['smallest_cluster']:
                        proto_galaxies = numpy.array(proto_cluster.galaxies)
                        ggroup = proto_galaxies[test]
                        include = True

                elif outlier_metric == 'vdisp':
                    proto_cluster = cluster(ggroup,-1)
                    vd = proto_cluster.vel_disp
                    vel = proto_cluster.velocity
                    vels = numpy.array([g.velocity for g in proto_cluster.galaxies])
                    dvels = numpy.array([abs(g.velocity-proto_cluster.velocity) for g in proto_cluster.galaxies])
                    test = numpy.array([dv < (nsigma*vd) for dv in dvels])
                    if countOf(test,True) >=self.selection['detection']['smallest_cluster']:
                        proto_galaxies = numpy.array(proto_cluster.galaxies)
                        ggroup = proto_galaxies[test]
                        include = True

                if plot:
                    from detections import plot_galaxies
                    plot_galaxies(proto_cluster.galaxies,colour='red')
                    xlims = plt.xlim()
                    ylims = plt.ylim()
                    accepted = numpy.array(proto_cluster.galaxies)
                    accepted = accepted[test]
                    plot_galaxies(accepted,colour='green',passback=True)
                    plt.xlim(xlims)
                    plt.ylim(ylims)
                    plt.title('%s %s'%(resample,outlier_metric))

                accepted = None
                if include:
                    accepted = cluster(ggroup,-1)

                    if verbose:
                        ot,ngb,nga,e0,e1 = outlier_metric,len(proto_cluster.galaxies),len(ggroup),degrees(proto_cluster.extent)*60.0,degrees(accepted.extent)*60.0,
                        print '%s %f-sigma cut reduced membership from %d to %d, and extent from %f to %f'%(ot,nsigma,ngb,nga,e0,e1)
                else:
                    if verbose:
                        ot,ngb,nga = outlier_metric,len(proto_cluster.galaxies),countOf(test,True)
                        print '%s %f-sigma cut reduced membership from %d to %d'%(ot,nsigma,ngb,nga)
                        

                if no_clobber and accepted == None:
                    accepted = cluster(ggroup,-1)
                    

                return accepted



    def perfect(self,data,halomass=13.5,verbose=True,nsigma=3.0,resample='bootstrap',outlier_metric='radius'):
        import numpy
        from tools import bootstrap as resampler
        
        if self.selection['mock_data'].has_key('halomass'):
            halomass = self.selection['mock_data']['halomass']

        if outlier_metric == False:
            outlier_metric = 'None'

        if resample == 'jacknife':
            from tools import jacknife as resampler
        try:
            dict = data
            kk = dict.keys()
            print 'Using data supplied as arg'
#            dict = self.cuts

        except:
            print 'Attempting to retrieve data'
            from dataIO import getData
            try:
                data = self.source_data
                ll = len(data['spatial']['x'])
            except:
                self.source_data = None
                dataset = self.selection['io']['dataset']
                catalogue = self.selection['io']['sourcefile']
                try:
                    obs = self.selection['io']['observed_data']
                except:
                    obs = False
                                        
                data = getData(catalogue,obs=obs,dataset=dataset)
                
                try:
                    self.source_data = data
                except:
                    pass
            self.get_input(data=self.source_data)
            dict = self.cuts

        galaxies = numpy.array(objectify(dict))
        source_population = galaxies       
        clusters = []
        names = []
        n = 0
        gals_names = {}
        
        hm = numpy.array([g.halomass for g in galaxies])
        filter = (hm >= halomass)
        galaxies = galaxies[filter]

        ids = unique([g.halo_id for g in galaxies])
        gals_names = {}
        clusters = {}
        for id in ids:
            clusters[id] = []
            
        for g in galaxies:
            g.vertices = []
            clusters[g.halo_id].append(g)

        perfect_clusters = []
        n = 0
        self.perfect_cluster_galaxies = []
        for candidate in clusters.keys():
            ggroup = clusters[candidate]
            if len(ggroup)>=self.selection['detection']['smallest_cluster'] and outlier_metric != 'None':
                out_test = self.perfect_outliers(ggroup,outlier_metric=outlier_metric,resample=resample,nsigma=nsigma)
                if out_test != None:
                    out_test.id = n
                    perfect_clusters.append(out_test)
                    self.perfect_cluster_galaxies = self.perfect_cluster_galaxies + list(ggroup)

#                    print ggroup[0].halo_id, perfect_clusters[-1].halo_id
#                    raise SystemExit(0)



                    n = n + 1
            elif len(ggroup)>=self.selection['detection']['smallest_cluster'] and outlier_metric == 'None':
                perfect_clusters.append(cluster(ggroup,n))
                self.perfect_cluster_galaxies = self.perfect_cluster_galaxies + list(ggroup)
                n = n + 1


        self.perfect_clusters = perfect_clusters
        self.perfect_cluster_galaxies = numpy.array(self.perfect_cluster_galaxies)
        self.perfect_field_galaxies = []
        
        source_population = numpy.array(source_population)
        hc = numpy.array([g.host_cluster == None for g in source_population])
        if len(hc) > 0:
            self.perfect_field_galaxies = source_population[hc]
        else:
            self.perfect_field_galaxies = []

        print '\t#Found %d perfect clusters, %d with N>=%d galaxies'%(len(clusters),len(perfect_clusters),self.selection['detection']['smallest_cluster'])
            

        return perfect_clusters
    


    def purity(self):
        purities = []
        nclus = len(self.clusters)
        if nclus == 0: return 0
        for cluster in self.clusters:
            ngal = len(cluster.galaxies)
            pure = [gals.false_detection for gals in cluster.galaxies]
            x = float(countOf(pure,True))
            y = float(len(pure))
            P = 1.0-(x/y)
            cluster.purity = P
            purities.append(P)
        self.pure = median(purities)
        return purities



    def quality(self):
        if len(self.clusters) == 0:
            return 0
        source_galaxies = objectify(self.cuts)
        qualities = []
        nclus = len(self.clusters)
        all_names = []

        for cluster in self.clusters:
            Ng = countOf([g.halo_id == cluster.halo_id for g in cluster.galaxies],True)
            Ngih = countOf([g.halo_id == cluster.halo_id for g in source_galaxies],True)
            Nfd = countOf([g.false_detection == True for g in cluster.galaxies],True)
            Ngic = len(cluster.galaxies)

            Q = (float(Ng)/float(Ngih))-(float(Nfd)/float(Ngic))
            
            qualities.append(Q)
            cluster.quality = Q
        self.qual = median(qualities)
        return qualities


    def quality_old(self):
        qualities = []
        nclus = len(self.clusters)
        if nclus == 0: return 0
        all_names = []
        for cluster in self.clusters:
            all_names = all_names + [g.uid for g in cluster.galaxies]

        for cluster in self.clusters:
            pure = [gals.false_detection for gals in cluster.galaxies]
            ngood = float(countOf(pure,False))
            nbad = float(countOf(pure,True))
            nspawn = countOf(all_names,cluster.name)
            Q = float(ngood-nbad)/float(nspawn)
            
            qualities.append(max(Q,0.0))
            cluster.quality = max(Q,0.0)
        self.qual = median(qualities)
        return qualities
  
    def accuracies(self):
        r = self.accuracy_r(self)
        s = self.accuracy_sm(self)
        v = self.accuracy_v(self)
        return [r,s,v]

    def accuracy_r(self):
        if len(self.perfect_clusters) == 0 and len(self.clusters) == 0: return 1.0
        if len(self.clusters) == 0: return 0.0
        perfect_names = [cluster.uid for cluster in self.perfect_clusters]
        accuracies = []
        mia = []
        nclus = len(self.clusters)
        for cluster in self.clusters:
            gal_radius = [hypot(g.x-cluster.x,g.y-cluster.y) for g in cluster.galaxies]
            gal_radius.sort()
            gal_radius = gal_radius[0:int(0.8*len(gal_radius))-1]
            r_det = mean(gal_radius)
            if cluster.name in perfect_names:
                i = indexOf(perfect_names,cluster.name)
                xs = [g.x for g in self.perfect_clusters[i].galaxies]
                ys = [g.y for g in self.perfect_clusters[i].galaxies]
                rs = [hypot( xs[j]-self.perfect_clusters[i].x, ys[j]-self.perfect_clusters[i].y) for j in xrange(len(xs))]
                r_per = mean(rs)
                a_r = 1.0-(abs(r_det-r_per)/(r_per+r_det))
                accuracies.append(a_r)
                cluster.acc_r = a_r
            else:
                if not cluster.name in self.mia:
                    self.mia.append(cluster.name)
        if len(accuracies) > 0: 
            self.acc_r = median(accuracies)
        else:
            accuracies.append(0.0)
            self.acc_r = 0.0
        return accuracies

    def accuracy_sm(self):
        if len(self.perfect_clusters) == 0 and len(self.clusters) == 0: return 1.0
        if len(self.clusters) == 0: return 0.0
        perfect_names = [cluster.uid for cluster in self.perfect_clusters]
        detected_names = [cluster.name for cluster in self.clusters]
        accuracies = []
        mia = []
        nclus = len(self.clusters)
        for cluster in self.clusters:
            sm_det = 0
            #make sure cluster isn't a substructure:
            if not ('sub' in cluster.uid):
                indices = compress([c.name == cluster.name for c in self.clusters],range(len(self.clusters)))
                sm_det = [self.clusters[i].stellarmass for i in indices]
                sm_det = sum(unique(sm_det))/(len(unique(sm_det))-countOf(unique(sm_det),0.0))
                if cluster.name in perfect_names:
                    i = indexOf(perfect_names,cluster.name)
                    sm_per = self.perfect_clusters[i].stellarmass
                    a_sm = 1.0-(abs(sm_det-sm_per)/(sm_per+sm_det))
                    accuracies.append(max(a_sm,0.))
                    cluster.acc_sm = max(a_sm,0.)
                else:
                    if not cluster.name in self.mia:
                        self.mia.append(cluster.name)
        if len(accuracies) > 0: self.acc_sm = median(accuracies)
        else:
            self.acc_sm = 0.0
            accuracies.append(0.0)
        return accuracies


    def accuracy_v(self):
        if len(self.perfect_clusters) == 0 and len(self.clusters) == 0: return 1.0
        if len(self.clusters) == 0: return 0.0
        perfect_names = [cluster.uid for cluster in self.perfect_clusters]
        detected_names = [cluster.name for cluster in self.clusters]
        accuracies = []
        mia = []
        nclus = len(self.clusters)
        for cluster in self.clusters:
            v_det = 0
            #make sure cluster isn't a substructure:
            if not ('sub' in cluster.uid):
                indices = compress([c.name == cluster.name for c in self.clusters],range(len(self.clusters)))
                v_det = []
                for i in indices:
#                    v_det = v_det + [abs(dv(gal.zcos,cluster.zcos)) for gal in self.clusters[i].galaxies]
                    v_det = v_det + [((gal.vx*gal.vx)+(gal.vy*gal.vy)+(gal.vz*gal.vz)) for gal in self.clusters[i].galaxies]
                v_det = sum(unique(v_det))/(len(unique(v_det))-countOf(unique(v_det),0.0))
                v_det = v_det**0.5
                if cluster.name in perfect_names:
                    i = indexOf(perfect_names,cluster.name)
#                    v_per = [abs(dv(gal.zcos,self.perfect_clusters[i].zcos)) for gal in self.perfect_clusters[i].galaxies]
                    v_per = [((gal.vx*gal.vx)+(gal.vy*gal.vy)+(gal.vz*gal.vz)) for gal in self.perfect_clusters[i].galaxies]
                    v_per = sum(unique(v_per))/(len(unique(v_per))-countOf(unique(v_per),0.0))
                    v_per = v_per**0.5
                    a_v = 1.0-(abs(v_det-v_per)/(v_per+v_det))
                    accuracies.append(max(a_v,0.))
                    cluster.acc_v = max(a_v,0.)
                else:
                    if not cluster.name in self.mia:
                        self.mia.append(cluster.name)
        if len(accuracies) > 0: self.acc_v = median(accuracies)
        else:
            self.acc_v = 0.0
            accuracies.append(0.0)
        return accuracies

    def fragmentation(self):
        if len(self.clusters) == 0: return 0.0
        all_names = [cluster.name for cluster in self.clusters]
        fragmentation = []
        for cluster in self.clusters:
            if not ('sub' in cluster.uid):
                indices = compress([c.name == cluster.name for c in self.clusters],range(len(self.clusters)))
                if len(indices) == 1:
                    f = 1.0
                else:
                    f = 1.0/float(len(indices))
                fragmentation.append(f)

        self.frag = median(fragmentation)
        return fragmentation

    def terminator(self,skynet=False):
        """
        How did the cluster detections generally terminate? Via a probthresh limit, or rho_crit?

        """

        if skynet: print '\tWe will need Sarah Connor for this method: going back in time...'
        terminator = 0
        terms = []
        for cluster in self.clusters:
            terms.append(cluster.truncation_history)

        if countOf(terms,'rho_crit') > countOf(terms,'probthresh'): terminator = 1
        if countOf(terms,'rho_crit') < countOf(terms,'probthresh'): terminator = -1

        return terminator
        

    def P_s(self,cluster_gals=None,field_gals=None,overdense_galaxies=None,plot=False):
        from tools import split_population,objectify,probthresh_filter,split_clusters
        from vperc import percolate,compare
        from parallel_qhull import qhull
        from parallel_vertices import deploy_vertices
        from parallel_detect import parallel_detect
        selection = self.selection
        if self.Ps != None:
            return self.Ps

        results = qhull(self.cuts)
        galaxies = deploy_vertices(results,self.cuts)[0]

        if not self.selection['detection'].has_key('mean_density'):
            selection['detection']['mean_density'] = parallel_detect(self.cuts,selection,set_density=True)

        if (cluster_gals == None) or (field_gals == None):
            cluster_gals,field_gals = split_population(galaxies,selection)

        if overdense_galaxies == None:
            overdense_galaxies = probthresh_filter(galaxies,selection['detection']['probthresh'],selection['detection']['mean_density'],lumweight=True)

        if not selection['mock_data'].has_key('fd0'):
            clusters = percolate(overdense_galaxies,selection['detection']['smallest_cluster'],selection['detection']['mean_density'],thresh=0.0,exclusions=False,real=False)
            tc,fd = split_clusters(clusters,selection,galaxies,pc_gals=cluster_gals,field_gals=field_gals)
            selection['mock_data']['fd0'] = len(fd)
            selection['mock_data']['tc0'] = len(tc)
        
        overdense_galaxies.sort(compare)
        clusters = percolate(overdense_galaxies,selection['detection']['smallest_cluster'],selection['detection']['mean_density'],thresh=selection['detection']['rho_crit'],exclusions=False,real=False)
        
        if self.clusters == None:
            self.addclusters(clusters)
        tc,fd = split_clusters(clusters,selection,galaxies,pc_gals=cluster_gals,field_gals=field_gals)
        selection['mock_data']['fd'] = len(fd)
        selection['mock_data']['tc'] = len(tc)
        try:
            self.Ps = 1.0-(float(selection['mock_data']['fd'])/float(selection['mock_data']['fd0']))
        except:
            self.Ps = 1.0
        return self.Ps
                          
    def ps1_mstar(self,halo=False,assoc=False,lf_correct=False,fix_inc=False,inc_mode='purity',blim=13.3,orca_inc='lum'):
        import numpy
        from dataIO import load
        for cl in self.clusters:
            self.stellarmass = numpy.nan
            self.esm = numpy.nan
            self.halomass = numpy.nan
            
        mstar_fit = load('/gal/r3/dmurphy/projects/panstarrs/ymstar.dat')
        lf_fit = None
        if lf_correct:
            lf_fit = load('/gal/r3/dmurphy/projects/panstarrs/ylfevoV3.dat')                    ### normal : evolving fit parameters with each redshift
                                                                                              ### V2: as V1, but resets fit to z=0 bJ LF
                                                                                              ### V3: as V1, but includes phi* evolution

            
        inc_data = None
        if fix_inc:
            inc_data = load('/gal/r3/dmurphy/projects/panstarrs/orca_%s.dat'%(inc_mode))

        inc = []
        if orca_inc == 'M*':
            print 'Applying ORCA incompleteness to stellar mass estimate (not luminosity)'
            inc = numpy.array([cl.get_luminosity(band='y',assoc=assoc,lf_data=lf_fit,lf_correct=lf_correct,fix_inc=fix_inc,inc_data=inc_data,blim=blim,get_inc=True) for cl in self.clusters])
            print 'Median / mean correction = %f,%f'%(median(inc),mean(inc))

        ylums = numpy.array([cl.get_luminosity(band='y',assoc=assoc,lf_data=lf_fit,lf_correct=lf_correct,fix_inc=fix_inc,inc_data=inc_data,blim=blim,apply_inc=len(inc)==0) for cl in self.clusters])

        zbins = mstar_fit.keys()
        zbins.sort()
        dz = zbins[1]-zbins[0]
        inbin = []
        for cl in self.clusters:
            try:
                inbin.append(int(cl.z/(zbins[1]-zbins[0])))
            except:
                inbin.append(numpy.nan)
#        inbin = [int(z/(zbins[1]-zbins[0])) for z in [cl.z for cl in self.clusters]]
        model = lambda P,x: (P[0]*x) + P[1]
        measured = []
        for m in xrange(len(self.clusters)):
            bin = inbin[m]
            ylum = ylums[m]
            if ylum > 0. and bin < len(mstar_fit.keys())-1:
                fit = mstar_fit[zbins[bin]]
                P = [fit['m'],fit['c']]
                mstar = model(P,ylum)
                if len(inc) > 0:
                    mstar = log10((10**mstar)*inc[m])


                self.clusters[m].stellarmass = mstar
                dl = [ylum-Le for Le in fit['error']['lum']]
                dla = [abs(ylum-Le) for Le in fit['error']['lum']]
                j = indexOf(dla,min(dla))
                if dl[j] != dl[-1]:
                    _dl = fit['error']['lum'][j+1]-fit['error']['lum'][j]
                    err = fit['error']['e'][j]*(1-(abs(ylum-fit['error']['lum'][j])/_dl))+fit['error']['e'][j+1]*(1-(abs(ylum-fit['error']['lum'][j+1])/_dl))
                else:
                    err = fit['error']['e'][j]
                self.clusters[m].esm = err
                measured.append(self.clusters[m])
            else:
                print ylum, bin

        print 'Have stellar mass estimates for %d of %d clusters'%(len(measured),len(self.clusters))
        hmeasured = []
        if halo:
#            halo_fit = load('/gal/r3/dmurphy/projects/panstarrs/halo_fit.dat')
            P = [ 0.9161919 , -0.80954158]
            P = [ 1.08609617,  0.94697272]
            model = lambda P,x: (P[0]*x) + P[1]
            for cl in measured:
                x = cl.stellarmass
                if x > 0:
                    cl.halomass = model(P,x)
                    hmeasured.append(cl)
            print 'Have halo mass estimates for %d of %d clusters'%(len(hmeasured),len(self.clusters))            
                    
                
        return measured



    def mock_purity(self,dz=0.02,plot=False):
        from science import science
        sci = science()
        z = arange(0,1.00,dz)
        self.clusters = numpy.array(self.clusters)
        self.cluster_galaxies = numpy.array(self.cluster_galaxies)
        purity = {}
        zz = []
        for _z in z:
            purity[_z] = 0.0
            zfilter = sci.filter(self,'zcos',_z,_z+dz)
            print len(zfilter.clusters)
            add = False
            if len(zfilter.clusters) > 0:
                purity[_z] = mean([c.hhf for c in zfilter.clusters])
                add = True
            if add:
                zz.append(_z)

        P = [purity[_z] for _z in zz]
        if plot:
            plt.plot(zz,P,'ko')
            plt.plot(zz,P,'k-')
        return zz,P
        
        
        

    def mock_completeness(self,dz=0.02,dhm=0.5,plot=False,sampling=5,mode='differential',counts=False,log=False,add_points=False,add_numbers=False,metric='completeness',maxBCG_mode=False,simple=False,cmap='default',count_thresh=1,smooth=False,norm=True,baseline=None,box_text_size=9):
        """\
        example usage:
        mock_completeness(plot=True,dhm=0.2,sampling=3,dz=0.05,counts=False,add_points=False,metric='sm',add_numbers=True,maxBCG_mode=False,simple=True,count_thresh=5,smooth=True)
        """

        if os.getenv('HOSTNAME') != 'dubris':
            import matplotlib.pyplot as plt
            import matplotlib
            if cmap == 'default':
                cmap=plt.cm.jet
            

        from science import science
        from dataIO import Progress
        import numpy
        from clusters import cluster
        from math import log10
        import warnings
        warnings.filterwarnings('ignore',module='numpy')
        if counts == True:
            counts = 'perfect'

        if metric == 'hhf':
            smooth = True

        if counts != False:
            smooth = False
            norm = False
        self.clusters = numpy.array(self.clusters)
        self.cluster_galaxies = numpy.array(self.cluster_galaxies)
        
        lookup = {}
        smooth_lookup = {}
        sci = science()
        z = arange(0,1.00,dz)
        h = arange(self.selection['mock_data']['halomass'],15.1,dhm)
        detected_haloes = unique([g.halo_id for g in self.cluster_galaxies])
        cluster_assigned_halos = unique([g.halo_id for g in self.clusters])

        detected_haloes = numpy.array([g.halo_id for g in self.cluster_galaxies])
        detected_stellar_mass = numpy.array([g.stellarmass for g in self.cluster_galaxies])

        cluster_assigned_halos = numpy.array([g.halo_id for g in self.clusters])

        perfect_detected_haloes = unique([g.halo_id for g in self.perfect_clusters])
        detected_dna = numpy.array([hash('%1.8f %1.8f'%(float(g.x),float(g.y))) for g in self.cluster_galaxies])

        perfect_dna = numpy.array([hash('%1.8f %1.8f'%(float(g.x),float(g.y))) for g in self.perfect_cluster_galaxies])
        stellar_mass = numpy.array([g.stellarmass for g in self.perfect_cluster_galaxies])
        pc_gal_hids = numpy.array([g.halo_id for g in self.perfect_cluster_galaxies])

        for c in self.perfect_clusters:
            c.detected = False

        performance = {}
        counter = {}
        
        xdata = []
        ydata = []
        zdata = []
        
        runs = len(z)*len(h)
        status = Progress(len(z))
        thresh = self.selection['detection']['smallest_cluster']
        zlist = []
        gridx = [_z for _z in z]
        gridy = [_h for _h in h]
        zlist = numpy.ones(runs,dtype='float')
        _h2 = h[5]
        _z2 = z[5]



        for _z in z:
            counter[_z] = {}
            for _h in h:
                counter[_z][_h] = 0

#        zlist = numpy.array([numpy.nan for i in arange(runs)])

        sequence = []
        for _z in z:
            performance[_z] = {}
            if mode == 'differential':
                zfilter = sci.filter(self,'z',_z,_z+dz)
            else:
                zfilter = sci.filter(self,'z',0,_z+dz)
            lookup[_z] = {}
            if len(zfilter.perfect_clusters) == 0:
                for _h in h:
#                    performance[_z][_h] = -0.5
                    sequence.append(0)
                    zlist[len(sequence)-1] = -0.5
                    lookup[_z][_h] = len(sequence)-1

            if len(zfilter.perfect_clusters) > 0:
                for _h in h:
                    sequence.append(0)
                    lookup[_z][_h] = len(sequence)-1
                    sm_det = 0.0
                    sm0 = 0.0
                    performance[_z][_h] = []                    
                    if mode == 'differential':
                        hfilter = sci.filter(zfilter,'halomass',_h,_h+dhm)
                    else:
                        hfilter = sci.filter(zfilter,'halomass',0,_h+dhm)
                    tot_sm = 0.0
                    det_sm = 0.0

                    
#                    condition = len(hfilter.perfect_clusters) >= count_thresh
                    condition = len(hfilter.perfect_clusters) >= 1

                    if condition:
                        halos_in_bin = [g.halo_id for g in hfilter.perfect_clusters]
                        if 1:
                            halos_in_bin = unique([g.halo_id for g in hfilter.perfect_clusters])
                            if metric == 'completeness':
                                for halo in halos_in_bin:
                                    if maxBCG_mode:
                                        thresh = 0.3
                                        pc = self.perfect_clusters[indexOf([g.halo_id for g in self.perfect_clusters],halo)]
                                        cluster_indices = numpy.where(detected_haloes == halo)[0]
                                        found = self.cluster_galaxies[cluster_indices]
                                        found = unique([f.parent for f in found])
                                        test = []
                                        for c in found:
                                            test.append(float(countOf([g.halo_id for g in c.galaxies],halo))/float(len(pc.galaxies)) >= 0.3)
                                            
                                        if True in test:
                                            performance[_z][_h].append(halo)

                                    else:
                                        mode2 = 'liberal'
#                                        mode2 = 'conservative'
                                        nn = {'conservative':countOf(cluster_assigned_halos,halo),'liberal':countOf(detected_haloes,halo)}
                                        n = nn[mode2]
                                        if mode2 == 'conservative' and n > 0:
                                            filter = numpy.array([c.halo_id == halo for c in self.clusters])
                                            clusters = self.clusters[filter]
                                            galaxies = []
                                            for c in clusters:
                                                galaxies = galaxies + c.galaxies
                                            g_hid = [g.halo_id for g in galaxies]
                                            n = countOf(g_hid,halo)

                                        if n >= thresh:
                                            performance[_z][_h].append(halo)
                                            pc = self.perfect_clusters[indexOf([g.halo_id for g in self.perfect_clusters],halo)]
                                            pc.detected = True
                                            
                                            


                            elif metric == 'sm':
                                if maxBCG_mode:
                                    filter = numpy.array([c.halo_id in halos_in_bin for c in self.clusters])
                                    matched = self.clusters[filter]
                                    sm = [c.stellarmass for c in matched]
                                    if len(sm) > 0:
                                        det_sm = sum(sm)
                                    else:
                                        det_sm = 0.0
                                    tot_sm = sum([c.stellarmass for c in hfilter.perfect_clusters])
                                    
                                    det_sm = []
                                    tot_sm = []
                                    for pc in hfilter.perfect_clusters:
                                        filter = numpy.array([c.halo_id == pc.halo_id for c in self.clusters])
                                        matched = self.clusters[filter]
                                        sm = [c.stellarmass for c in matched]
                                        if len(sm) > 0:
                                            det_sm.append(sum(sm))
                                        else:
                                            det_sm.append(0.0)
                                        tot_sm.append(pc.stellarmass)

                                else:
                                    pcs = {}
                                    for halid in halos_in_bin:
                                        pcs[halid] = []
                                    for c in hfilter.perfect_clusters:
                                        pcs[c.halo_id].append(c)
                                    new_cls = []
                                    for c in pcs.keys():
                                        gals = []
                                        for sub_det in pcs[c]:
                                            gals = gals + list(sub_det.galaxies)
                                        cl = cluster(gals,0)
                                        new_cls.append(cl)
                                    for c in new_cls:
                                        tot_sm = tot_sm + c.stellarmass
                                        cdna = numpy.array([hash('%1.8f %1.8f'%(float(g.x),float(g.y))) for g in c.galaxies])
                                        for dna in cdna:
                                            if dna in detected_dna:
                                                det_sm = det_sm + c.galaxies[indexOf(cdna,dna)].stellarmass

                    if metric == 'hhf':
#                    if len(hfilter.clusters) > 0 and metric == 'hhf':
                        if 1:
                            if len(hfilter.clusters) > 0 and len(hfilter.perfect_clusters) > 0:
                                performance[_z][_h] = [_c.hhf for _c in hfilter.clusters]
                            elif len(hfilter.perfect_clusters) > 0:
                                performance[_z][_h] = [0.0]


#                            else:
#                                performance[_z][_h] = [3.0/5.0]

                    if condition or metric == 'hhf':
#                    if condition or (len(hfilter.clusters) > 0 and metric == 'hhf'):
                        if metric == 'completeness':
                            performance[_z][_h] = float(len(performance[_z][_h]))/float(len(halos_in_bin))
                        elif metric == 'sm':
                            if maxBCG_mode:
                                if len(det_sm) > 0:
                                    sample = []
                                    for i in xrange(len(det_sm)):
                                        try:
                                            p = (log10(float(det_sm[i]))-log10(float(tot_sm[i])))+1
                                        except ValueError:
                                            p = -1.0
                                        sample.append(p)
                                    performance[_z][_h] = median(sample)

                                    performance[_z][_h] = 1-(abs(1-performance[_z][_h]))
                                    
                                    if performance[_z][_h] == -1:
                                        performance[_z][_h] = 0.0

                                else:
                                    try:
                                        performance[_z][_h] = (log10(float(det_sm))-log10(float(tot_sm)))+1
                                    except:
                                        performance[_z][_h] = -1.0

                            else:
                                performance[_z][_h] = float(det_sm)/(tot_sm)

                        elif metric == 'hhf' and len(performance[_z][_h]) > 0:
#                            if len(hfilter.clusters) > 0 and len(hfilter.perfect_clusters) > 0:
                            if len(hfilter.perfect_clusters) > 0:
#                                performance[_z][_h] = mean(performance[_z][_h])
                                performance[_z][_h] = median(performance[_z][_h])
                            else:
                                performance[_z][_h] = -0.5
                        elif metric == 'hhf' and len(performance[_z][_h]) == 0:
                            performance[_z][_h] = -0.5
                                
                            
                        if counts == False:
                            if type(performance[_z][_h]) != type([]):
                                zdata.append(performance[_z][_h])
                                zlist[len(sequence)-1] = performance[_z][_h]
                                xdata.append(_z)
                                ydata.append(_h)
                                counter[_z][_h] = len(hfilter.perfect_clusters)
                        
                        elif counts == 'perfect':
                            zdata.append(float(len(hfilter.perfect_clusters)))
                            zlist[len(sequence)-1] = float(len(hfilter.perfect_clusters))
                            xdata.append(_z)
                            ydata.append(_h)

                        elif counts == 'detected':
                            zdata.append(float(len(hfilter.clusters)))
                            zlist[len(sequence)-1] = float(len(hfilter.clusters))
                            xdata.append(_z)
                            ydata.append(_h)
                    
                    else:
                        performance[_z][_h] = -0.5
                        zlist[len(sequence)-1] = -0.5
            status.update(0)


        if smooth:
            smoothed = {}
            seq = []
            zl = numpy.ones(runs,dtype='float')
            _h2 = h[5]
            _z2 = z[5]


            for _z in z:
                smoothed[_z] = {}
                for _h in h:
                    seq.append(0)
                    try:
                        smooth_lookup[_z][_h] = len(seq)-1
                    except:
                        smooth_lookup[_z] = {}
                        smooth_lookup[_z][_h] = len(seq)-1                        

                    
#                    smoothed[_z][_h] = -0.5
                    zl[len(seq)-1] = -0.5
                    ref = arange(-1,1.1,1,dtype='int')
                    i = indexOf(z,_z)
                    j = indexOf(h,_h)
                    n_grid = {}
                    dat_grid = {}
                    for r1 in ref:
                        for r2 in ref:
                            try:
                                if performance[z[r1+i]][h[r2+j]] >= 0 and performance[z[i]][h[j]] >= 0:
                                    try:
                                        n_grid[(r1,r2)] = counter[z[r1+i]][h[r2+j]]
                                        dat_grid[(r1,r2)] = performance[z[r1+i]][h[r2+j]]
                                    except:
                                        pass
                            except KeyError:
                                pass
                            except IndexError:
                                pass
                    
                    try:
                        if sum(n_grid.values()) >= count_thresh:
                            smoothed[_z][_h] = sum(dat_grid.values())/len(dat_grid.values())
#                        smoothed[_z][_h] = sum(dat_grid.values())/sum(n_grid.values())
                            zl[len(seq)-1] = smoothed[_z][_h]
                    except:
                        pass

                        
                        
            print 'Smoothing done'
            zlist = zl
            performance = smoothed
#            lookup = smooth_lookup
            

            min_thresh = 0.0
            if baseline != None:
                min_thresh = -99.0
                zlv2 = numpy.ones(runs,dtype='float')-1000.0
                new_performance = {}
                _matched = []
                _unmatched = []
                for zc in performance.keys():
                    if baseline.has_key(zc):
                        for hm in performance[zc].keys():
                            if baseline[zc].has_key(hm):
                                j = smooth_lookup[zc][hm]
#                                k = lookup[zc][hm]
                                zlv2[j] = zlist[j]-baseline[zc][hm]
                                _matched.append([zc,hm])
                                try:
                                    new_performance[zc][hm] = zlv2[j]
                                except:
                                    new_performance[zc] = {}
                                    new_performance[zc][hm] = zlv2[j]
                            else:
                                _unmatched.append([zc,hm])
                    else:
                        for hm in performance[zc].keys():
                            _unmatched.append([zc,hm])
                zlist = zlv2
                if len(_unmatched) > 0:
                    print 'There were %d cells unmatched between this analysis and the benchmark. They were:'%(len(_unmatched))
                    for u in _unmatched:
                        print '(z,MH)=(%1.2f,%2.2f)'%(u[0],u[1])
                performance = new_performance

        if plot:
            from matplotlib.mlab import griddata
            from numpy import linspace,log10
            if counts != False:
                metric = counts
            labels = {'completeness':'Completeness','sm':'M* detection fraction',\
                      'perfect':'Perfect cluster counts','detected':'Detected cluster counts',\
                      'hhf':'Purity'}



            if simple:
                import numpy
                _gridx = list(numpy.array(gridx) - (0.5*dz))
                _gridy = list(numpy.array(gridy) - (0.5*dhm))
                val2 = numpy.array(zlist).reshape(len(_gridx),len(_gridy))
                val2 = val2.transpose()
                masked_array = numpy.ma.array(val2, mask=val2<min_thresh)

                x_grid,y_grid = numpy.meshgrid(_gridx,_gridy)
#                vals = numpy.array(z)
#                val2 = vals.reshape(len(xdata),len(ydata))
#                ax = plt.pcolormesh(x_grid,y_grid,val2,cmap=plt.cm.binary)

                mynorm = None
                if norm == True:
                    if baseline != None:
                        zlist = numpy.array(zlist)
                        filter = (zlist > min_thresh)
                        _zl = zlist[filter]
                        minval = min(_zl)
                        maxval = max(_zl)
                        vv = max(abs(minval),abs(maxval))
                        mynorm = plt.Normalize(vmin=-1.0*vv,vmax=vv)
                        cmap = plt.cm.RdBu
#                        cmap = plt.cm.Spectral
                    else:
                        mynorm = plt.Normalize(vmin=0.0,vmax=1.0)
                    print 'Normalising to peak values'
                ax = plt.pcolormesh(x_grid,y_grid,masked_array,cmap=cmap,norm=mynorm)
#                cb = plt.colorbar(plt.gca(),format='%.1f')


            else:
                samplingx = len(z)*sampling
                xi = linspace(min(xdata),max(xdata),samplingx)
                yi = linspace(min(ydata),max(ydata),samplingx)
                if log and counts != False:
                    zdata = log10(zdata)
                zi = griddata(xdata,ydata,zdata,xi,yi)
                CS = plt.contourf(xi,yi,zi,250,cmap=plt.cm.jet)

            if add_points:
                plt.plot(xdata,ydata,'ko')

            if add_numbers and simple:
                for _x in gridx:
                    i = indexOf(gridx,_x)
                    for _y in gridy:
                        j = indexOf(gridy,_y)
                        try:
                            val = float(masked_array[j][i])
                        except:
                            continue
                        if val > 0.94:
                            txt_colour = 'w'
                        elif val > 0.30:
                            txt_colour = 'k'
                        else:
                            txt_colour = 'w'

                        if cmap in [plt.cm.RdBu,plt.cm.Spectral]:
                            if abs(val) < 0.20:
                                txt_colour = 'k'
                            else:
                                txt_colour = 'w'

                        if numpy.isnan(val):
                            continue
                        val = '%1.2f'%(val)
                        plt.text(_x,_y,val,fontsize=box_text_size,horizontalalignment='center',verticalalignment='center',color=txt_colour)


            elif add_numbers:
                for i in xrange(len(xdata)):
                    xpos,ypos = xdata[i],ydata[i]
                    if counts:
                        val = '%1.2f'%(zdata[i])
                    else:
                        try:
                            val = '%1.2f'%(performance[xpos][ypos])
                        except:
                            val = ''
                            pass
                    txt_col = 'white'
                    if zdata[i] <= 0.7:
                        txt_col = 'k'

                    if 1:
#                    if xpos < max(z) and xpos > min(z) and ypos < max(h) and ypos > min(h):                        
                        plt.text(xpos,ypos,val,fontsize=8,horizontalalignment='center',color=txt_col)

            plt.xlim(0,max(z))
            plt.ylim(min(h),max(h))
#            plt.ylim(min(ydata),max(ydata))
            plt.ylabel('Halo mass')
            plt.xlabel('Redshift (z)')
            cb = plt.colorbar(format='%.1f')
            title = labels[metric]
            if log:
                title = title+' [log10]'
            if maxBCG_mode:
                title = 'maxBCG-defined '+title
            cb.ax.set_ylabel(title)
                
        try:
            performance = new_performance
        except:
            pass
        return performance


    def mock_smfrac(self,dz=0.02,dhm=0.5,plot=False,sampling=5,mode='differential'):
        from science import science
        from dataIO import Progress
        import numpy
        sci = science()
        z = arange(0,1.00,dz)
        h = arange(self.selection['mock_data']['halomass'],15.1,dhm)
        detected_haloes = unique([g.halo_id for g in self.cluster_galaxies])
        detected_dna = [hash('%1.8f %1.8f'%(float(g.x),float(g.y))) for g in self.cluster_galaxies]
        completeness = {}
        stellarmass = {}
        xdata = []
        ydata = []
        zdata = []
        
        runs = len(z)*len(h)
        status = Progress(len(z))

        for _z in z:
            stellarmass[_z] = {}
            if mode == 'differential':
                zfilter = sci.filter(self,'z',_z,_z+dz)
            else:
                zfilter = sci.filter(self,'z',0,_z+dz)
            if len(zfilter.perfect_clusters) > 0:
                for _h in h:
                    tot_sm = 0.0
                    det_sm = 0.0
                    det_frac = []
                    stellarmass[_z][_h] = []
                    if mode == 'differential':
                        hfilter = sci.filter(zfilter,'halomass',_h,_h+dhm)
                    else:
                        hfilter = sci.filter(zfilter,'halomass',0,_h+dhm)
                    if len(hfilter.perfect_clusters) > 0:
                        #define detected?
                        for c in hfilter.perfect_clusters:
                            tot_sm = tot_sm + c.stellarmass
                            cdna = [hash('%1.8f %1.8f'%(float(g.x),float(g.y))) for g in c.galaxies]
                            n = 0
                            matches = []
                            det_obs = []
                            for dna in cdna:
                                if dna in detected_dna:
                                    matches.append(self.cluster_galaxies[indexOf(detected_dna,dna)].parent)
                                    det_obs.append(c.galaxies[indexOf(cdna,dna)].stellarmass)
                                    det_sm = det_sm + c.galaxies[indexOf(cdna,dna)].stellarmass
                                    n = n + 1
                            if len(matches) > 0:
                                matches = unique(matches)
                                det_frac.append(float(sum(det_obs))/float(c.stellarmass))
#                                det_sm = sum([c.stellarmass for c in matches])
                                
                        stellarmass[_z][_h] = float(det_sm)/float(tot_sm)
#                        stellarmass[_z][_h] = mean(det_frac)
                        xdata.append(_z)
                        ydata.append(_h)
                        zdata.append(stellarmass[_z][_h])
                    else:
                        stellarmass[_z][_h] = -0.5
            status.update(0)

        if plot:
            from matplotlib.mlab import griddata
            from numpy import linspace
#            ydata = [10**y for y in ydata]
            samplingx = len(z)*sampling
            xi = linspace(min(xdata),max(xdata),samplingx)
            yi = linspace(min(ydata),max(ydata),samplingx)
            zi = griddata(xdata,ydata,zdata,xi,yi)
            CS = plt.contourf(xi,yi,zi,250,cmap=plt.cm.jet)
            plt.xlim(0,max(z))
            plt.ylim(min(h),max(h))
#            plt.ylim(min(ydata),max(ydata))
            plt.ylabel('Halo mass')
            plt.xlabel('Redshift (z)')
            cb = plt.colorbar(format='%.1f')
            plt.ylabel(r'log$_{10}$M$_{halo}$ (M$_{\odot}$)')
#            cb.ax.set_ylabel(r'$M_{*}^{det}$/$M_{*}^{true}$')
            cb.ax.set_ylabel(r'M$_{*}$detection fraction')
        return stellarmass



    def completeness(self):
        if len(self.perfect_clusters) == 0 and len(self.clusters) == 0: return 1.0
        perfect_names = [pc.halo_id for pc in self.perfect_clusters]
        detected_names = unique([g.halo_id for g in self.cluster_galaxies])

        matches = [name in detected_names for name in perfect_names]
        matches = countOf(matches,True)
        try:
            self.comp = float(matches)/float(len(perfect_names))
        except ZeroDivisionError:
            self.comp = 0.0
        return self.comp
        

    def completeness_old(self,mlow=13.5,dm=0.1,rlimit=24):
        from tools import all_good_detcted, count_names
        if len(self.perfect_clusters) == 0 and len(self.clusters) == 0: return 1.0
        massbins = arange(mlow,15.5,dm)
        good_detected_names = all_good_detcted(self.clusters,mlow)

        redshifts,halomass = [],[]
        cullnames = []
        for cluster in self.perfect_clusters:
            mag_check = [galaxy.rmag < rlimit for galaxy in cluster.galaxies]
            if countOf(mag_check,True) >=5:
                cluster.what_is_bcg()
                redshifts.append(cluster.bcg.zcos)
                halomass.append(cluster.halomass)
                cullnames.append(cluster.uid)

        completeness = count_names(good_detected_names, cullnames)
        try:
            completeness = 1.0*completeness/len(cullnames)
        except ZeroDivisionError:
            #means we have no perfect clusters, but have detected clusters: defaults completeness to 0.0
            completeness = 0.0
        self.comp = completeness
        return completeness
    
    def addclusters(self,clusterlist,metrics=False,zdist=False,perfect=False,data=None,verbose=True,outlier_metric='radius',force_outlier=False,gcross=True):
        """\
        det.addclusters(self,clusterlist,metrics=False,zdist=False,perfect=False,data=None,verbose=True,outlier_metric='radius',force_outlier=False,gcross=True) => void
        \t Add ORCA-detected clusters to the detection class. This will build in the object-inheritance as well physically add them to the list.

        clusterlist     :  the list of cluster object
        metrics         :  if this is a mock, calculate the purity, completeness etc for the detection
        zdist           :  create redshift histogram data => det.redshift_dist
        perfect         :  these clusters are "perfect" clusters derived from the DM simulation data (not ORCA)
        data            :  source galaxy catalogue data (stored as a dict) => det.cuts
        verbose         :  display information about the added clusters
        outlier_metric  :  if ORCA must eject outlying members, what is the variable by which they are expelled?
        force_outlier   :  eject outlying galaxy members?
        gcross          :  cross-check all galaxies within detection to determine if they are shared with > 1 cluster


        """


        import numpy
        emptyarr = numpy.array([])

        if type(perfect) == type(emptyarr):
            perfect = list(perfect)

        if clusterlist == None:
            self.clusters = []
            return
        if len(clusterlist) == 2:
            if ((type(clusterlist[0]) == type([]) or type(clusterlist[0]) == type(emptyarr))):
                clusterlist,perfect = clusterlist
        cluster_galaxies = []
        self.clusters = clusterlist
        for i in xrange(len(self.clusters)):
            c = self.clusters[i]
            c.position = i
            c.parent = self
            if not c.sub_id:
                try:
                    c.sub_id = self.selection['detection']['subset_index']
                except:
                    pass
        

            if c.selection == None:
                c.photometric_selection(verbose=verbose)
        if len(self.clusters) > 0:
            if sum(self.clusters[0].merge_crits.values()) > 0:
                for key in self.merge_crits.keys():
                    self.merge_crits[key] = self.clusters[0].merge_crits[key]

            if self.clusters[0].survey_mean_density != None:
                md = self.clusters[0].survey_mean_density
                if verbose: print '\tUpdating mean_density %f -> %f'%(self.selection['detection']['mean_density'],md)

            for cluster in self.clusters:
                cluster.merge_crits = {1:0,2:0,3:0,4:0,5:0}
                for galaxy in cluster.galaxies:
                    cluster_galaxies.append(galaxy)
                    galaxy.parent = cluster
                    galaxy.position = indexOf(cluster.galaxies,galaxy)
        self.cluster_galaxies = cluster_galaxies
        if not cluster_galaxies:
            self.cluster_galaxies = []
        
        try:
            acg = []
            for c in self.clusters:
                acg = acg + list(c.associated_members)
            self.associated_cluster_galaxies = acg
        except:
            pass

        try:
            uber_clusters = []
            uber_cluster_galaxies = []
            for c in self.clusters:
                ucl = c.uber_cluster
                ucl.parent = self
                uber_clusters.append(ucl)
                for g in ucl.galaxies:
                    if (g in c.galaxies) == False:
                        uber_cluster_galaxies.append(g)
            self.uber_clusters = uber_clusters
            self.uber_cluster_galaxies = uber_cluster_galaxies
        except:
            pass
                        


        if 'mock' in self.selection['io']['dataset'] and perfect != False:
            if perfect == True:
                if not data:
                    from dataIO import getData
                    try:
                        data = self.cuts
                        ll = len(data['spatial']['x'])
                    except:
                        dataset = self.selection['io']['dataset']
                        catalogue = self.selection['io']['sourcefile']
                        obs = self.selection['io'].has_key('observed_data')
                        data = getData(catalogue,obs=obs,dataset=dataset)
                    try:
                        self.cuts = data
                    except:
                        pass

                    if self.selection['mock_data'].has_key('outlier_metric'):
                        outlier_metric = self.selection['mock_data']['outlier_metric']
                        print 'Overriding outlier metric to <<%s>>'%(outlier_metric)

                    perfect = self.perfect(data,self.selection['mock_data']['halomass'],outlier_metric=outlier_metric)
            else:
                if force_outlier:
                    print 'Applying outlier test to existing perfect cluster set'
                    if self.selection['mock_data'].has_key('outlier_metric'):
                        print 'Overriding outlier metric to <<%s>>'%(outlier_metric)
                    pc = []
                    for candidate in perfect:
                        out = self.perfect_outliers(candidate.galaxies,outlier_metric=outlier_metric,verbose=False)
                        if out != None:
                            pc.append(out)
                    print '\nRejected %d/%d perfect clusters'%(len(perfect)-len(pc),len(perfect))
                    perfect = pc

                    
                self.perfect_clusters = perfect
            perfect_cluster_galaxies = []

            for i in xrange(len(self.perfect_clusters)):
                c = self.perfect_clusters[i]
                c.position = i
                c.parent = self
                perfect_cluster_galaxies = perfect_cluster_galaxies + list(c.galaxies)
                if c.selection == None:
                    c.photometric_selection(verbose=verbose)
            self.perfect_cluster_galaxies = perfect_cluster_galaxies

        if zdist:
            self.redshift_dist = self.zdist()
        if(metrics) and len(self.clusters) > 0:
            self.purity()
            self.completeness()
            self.accuracy_r()
            self.accuracy_sm()
            self.accuracy_v()
            self.quality()
            self.fragmentation()
            self.terminator()
        try:
            sci = science()
            sci.richness(self)
        except NameError:
            try:
                import science
            except ImportError:
                print '\tCannot import science, not computing richness...'
                return
            sci = science.science()
            try:
                if len(self.clusters) > 0:
                    sci.richness(self)
                else:
                    for c in self.clusters:
                        c.richness = 1.0
                        c.exact_richness = 1.0


            except:
                for c in self.clusters:
                    c.richness = 1.0
                    c.exact_richness = 1.0

        

        if gcross:
            self.gcross()
        return
        

    def gcross(self,verbose=False):
        if verbose: print 'Performing galaxy cross-checks...'
        gids = {}
        gals = list(self.cluster_galaxies)
        try:
            assoc = list(self.associated_cluster_galaxies)
        except:
            assoc = []
        gals = gals + list(assoc)

        for g in gals:
            gids[g.objID] = {'member':[],'assoc':[]}

        for g in gals:
            cl_id = g.parent.dna
            if g in assoc:
                gids[g.objID]['assoc'].append(cl_id)
            else:
                gids[g.objID]['member'].append(cl_id)
                     

        for g in gals:
            g.memberOf = gids[g.objID]
        if verbose:print '....done'

    def metric_plots(self,plot=True,dz=0.05,zmax=1.5,mhalo=[],passback=False,variable='',legend_label='',style=''):
        Names = ['Purity','Completeness','Quality','Accuracy(r)','Accuracy(M*)',\
                     'Accuracy(v)','Fragmentation']
        if variable != '':
            if not variable in Names:
                for name in Names:
                    if variable in name:
                        variable = name
                        break
                if not variable in Names:
                    print '\tMetric %s not recognised - please re-submit using:\n %s'%(variable,str(Names))
                    return
        import science
        sci = science.science()
        zrange = arange(0.0,zmax+dz/2.0,dz)
        results = []
#        mhalo = [0.5e14,1e14,2e14,10e14,20e14]
        if len(mhalo) > 0:
            results = []
            cutresults = {}
            for m in mhalo:
                cutresults[m] = []

            for m in mhalo:
                mfilter = sci.filter(self,'halomass',10**log10(m),10**(log10(m)+10))
                for z in zrange:
                    filter = sci.filter(mfilter,'z',z,z+dz)
                    metric = sci.metrics(filter)
                    if variable != '':
                        thisz = {'z':z,variable:mean(metric[variable])}
                    else:
                        thisz = {'z':z+(dz/2.0),'C':metric['Completeness'],'P':mean(metric['Purity'])}
                        thisz['Q'] = mean(metric['Quality'])
                        thisz['A_v'] = mean(metric['Accuracy(v)'])
                    cutresults[m].append(thisz)
            results = cutresults
            #.append(cutresults[m])
            
        else:
            cutresults = []
            for z in zrange:
                filter = sci.filter(self,'z',z,z+dz)
                metric = sci.metrics(filter)
                if variable != '':
                    thisz = {'z':z,variable:mean(metric[variable])}
                else:
                    thisz = {'z':z,'C':metric['Completeness'],'P':mean(metric['Purity'])}
                    thisz['Q'] = mean(metric['Quality'])
                    thisz['A_v'] = mean(metric['Accuracy(v)'])
                cutresults.append(thisz)
            results.append(cutresults)

                
        if plot:
            if len(results) == 1:
                colours = ['b','g','r','c','m','y','k']
                colour = iter(colours)
                z = []
                P = []
                C = []
                Q = []
                Av = []
                generic = []
                for i in xrange(len(results[0])):
                    if variable != '':
                        z.append(results[0][i]['z'])
                        generic.append(results[0][i][variable])
                    else:
                        z.append(results[0][i]['z'])
                        P.append(results[0][i]['P'])
                        C.append(results[0][i]['C'])
                        Q.append(results[0][i]['Q'])
                        Av.append(results[0][i]['A_v'])

                if variable != '':
                    if legend_label != '':
                        plt.plot(z,generic,style,label=legend_label)
                    else:
                        plt.plot(z,generic,style)
                    if passback:
                        return
                    else:
                     plt.xlim(0,zmax)
                     plt.xlabel('z')
                     plt.ylabel('Performance')
                     plt.title(variable)
                     plt.show()   
                     return
                        
                
                plt.plot(z,C,'-k',color=colour.next(),label='Completeness')
                plt.plot(z,P,color=colour.next(),label='Purity')
                plt.plot(z,Q,color=colour.next(),label='Quality')
                plt.plot(z,Av,color=colour.next(),label='A_<v>')
                plt.ylim(0.0,1.0)
                plt.xlim(0.0,1.5)
                plt.text(0.95,0.70,self.print_sel(),fontsize=10,\
                         bbox=dict(facecolor='red', alpha=0.5),va='baseline')
                plt.xlabel('redshift')
                plt.ylabel('Performance')
                plt.legend()
                plt.show()
            else:
                colours = ['b','g','r','c','m','y','k']
                colour = iter(colours)
                for m in mhalo:
                    style = iter(['-','--','-.',':'])
                    col = colour.next()
                    z = []
                    P = []
                    C = []
                    Q = []
                    Av = []
                    for i in xrange(len(results[m])):
                        z.append(results[m][i]['z'])
                        P.append(results[m][i]['P'])
                        C.append(results[m][i]['C'])
                        Q.append(results[m][i]['Q'])
                        Av.append(results[m][i]['A_v'])
                
                    plt.plot(z,C,linestyle=style.next(),color=col,label='Completeness')
#                    plt.plot(z,P,linestyle=style.next(),color=col,label='Purity')
#                    plt.plot(z,Q,linestyle=style.next(),color=col,label='Quality')
#                    plt.plot(z,Av,linestyle=style.next(),color=col,label='A_<v>')

                plt.ylim(0.0,1.0)
                plt.xlim(0.0,1.5)
                plt.xlabel('redshift')
                plt.ylabel('Performance')
                plt.text(0.95,0.80,self.print_sel(),fontsize=10,\
                         bbox=dict(facecolor='red', alpha=0.5),va='baseline')

#                plt.legend()
                plt.show()
        return
                    


    def co_detect_colour(self,plot=True):
        """\
        detections.co_detect_colour(self,plot=True) -> void
    
        \t For a detections class derived from the co-detection of
        clusters in two colours - base (ie, the current colour) and next
        (ie, the colour a cluster must also have been detected in).

        Plots the dispersion of cluster colours in the next colour
        compared to the base colour.

        Eg. For a detection in g-r, will plot the g-r colour for each
        detection slice vs the mean and width of the cluster colour
        distribution in r-i.     

        """

        if not self.colour_colour:
            print '\tNo co-detections made with this detection instance'
            return None
        cc = self.colour_colour
        if len(cc['this']['top']) == 1: plot=False
        
        x = cc['this']
        y = cc['next']
        if not plot: return [x,y]
        #best plot is mean of this colour vs spread in next colour, so...
        if plot:
            data = [x,y]
            x = data[0]['middle']
            y1 = data[1]['bottom']
            y2 = data[1]['middle']
            y3 = data[1]['top']
            error = []
            error.append(y1)
            error.append(y3)
            plt.errorbar(x,y2,yerr=error,marker='o',color='b',linestyle='None')

            plt.xlabel(cc['this']['name'])
            plt.ylabel(cc['next']['name'])
            plt.show()
        return

    def magwidth(self,colour='default',plot=False,dmag=0.2,poly=2):
        """\
        detections.magwidth(self,colour='default',plot=False,dmag=0.2,poly=2) -> numpy.poly1d

        Returns the width function in arg:colour for a detections class, using data
        binned in magnitude widths arg:dmag, fitted to a polynomial of order arg:poly.

        #Process:
        #1. Rotate to remove CMR slope (fit from either earlier fit or on-the-fly)
        #2. For each cluster, mean intercept -> translate data to c=0
        #3. Asymmeterise: mirror all c<0 points
        #4. stack all cluster cmr's
        #5. Bin the galaxies in magnitude, with width dmag
        #6. For each bin, find maximum colour (or av/median of N highest/other prescription)
        #7. Fit over magnitude : linear, N-even polynomial, spline, exponential (currently only polynomial)
        #8. By-eye follow-up or goodness of fit in SciPy..... (if required)
        """
        if colour == 'default':
            try:
                colour = self.selection['defaults']['colour']
            except:
                colour = 'gr'


        import math
        import numpy

        # need to find out if there was a gradient set to the previous colour...            
        rotate = False
        colour_list = ['gr','ri','iz','zy']
        magnitude = colour[-1]+'mag'
        limits = {'ri':{'xlim':[15.5,24],'ylim':[-0.2,2.0]},'gr':{'xlim':[15.5,24],'ylim':[-0.2,2.3]},\
                  'iz':{'xlim':[15.5,24],'ylim':[-0.2,1.2]},'zy':{'xlim':[15.5,24],'ylim':[-0.2,0.6]}}
        if colour == colour_list[0]:
            prev_colour = colour
            if self.selection['colours'][colour]['slope'] != 0: 
                rotate = True
                print '\tUsing this (%s) colour gradient to rotate data'%(colour)

        else:
            prev_colour = ''
            for i in xrange(len(colour_list)):
                if colour == colour_list[i]:
                    try:
                        prev_colour = colour_list[i-1]
                        if self.selection['colours'][prev_colour]['slope'] != 0: 
                            rotate = True
                            print '\tUsing previous (%s) colour gradient to rotate data'%(prev_colour)
                    except:
                        print '\tNot expecting anything to happen here'
                    #first colour, so ignore...


        if rotate and self.selection['colours'][prev_colour]['fitfuncs'].has_key('gradient_function'):
            print '\t#detections.cmr(): detected %s gradient function - fitting slope from current clusters'%(prev_colour)
            self.selection['colours'][prev_colour]['slope'] = average(self.cmr(fit=True,colour=prev_colour))


        if not rotate:
            #need to find the mean slope for this detection...
            gradient = average(self.cmr(fit=True,colour=prev_colour))
            print '\tFitting on the fly from mean CMR slope'
            
        else:
            gradient = self.selection['colours'][prev_colour]['slope']            


        stacked_m = []
        stacked_c = []

        clusters = copy(self.clusters)
        for cluster in clusters:
            if not rotate: gradient = cluster.cmr(fit=True)[0]
            allm = []
            allc = []
            col = 'galaxy.'+colour
            mag = 'galaxy.'+magnitude
            for galaxy in cluster.galaxies:
                allm.append(eval(mag))
                allc.append(eval(col))
        
            meanmag = mean(allm)
            meanc = mean(allc)
            y0 = 0.0
            x0 = (-1.0*meanc/gradient)+meanmag
            r = []
            angle = []
            x_old = []
            y_old = []
            bcgx_old = []
            bcgy_old = []
            bcg_r = []
            bcg_ang = []
            for i in xrange(len(allm)):
                allm[i] = allm[i] - x0
                allc[i] = allc[i] - y0
                r.append(((allm[i]**2)+(allc[i]**2))**0.5)
                x_old.append(allm[i]/r[i])
                y_old.append(allc[i]/r[i])
                angle.append(math.atan2(allc[i],allm[i]))

            fit_angle = math.atan2(-1.0*gradient,-1.0)
            #need to rotate by 2pi - fit_angle
            d_theta = math.pi - fit_angle
#            print '\tRotating by %f'%(d_theta)
            x_new = []
            y_new = []
            bcgx_new,bcgy_new = [],[]
            for i in xrange(len(x_old)):
                angle[i] = angle[i] + d_theta
                x_new.append(math.cos(angle[i]))
                y_new.append(math.sin(angle[i]))

            for i in xrange(len(x_old)):
                x_new[i] = r[i]*x_new[i]
                y_new[i] = r[i]*y_new[i]
                x_old[i] = r[i]*x_old[i]
                y_old[i] = r[i]*y_old[i]
                
                x_new[i] = x_new[i] + x0 
                y_new[i] = y_new[i] + y0 + meanc
                x_old[i] = x_old[i] + x0
                y_old[i] = y_old[i] + y0

#            plt.plot(x_old,y_old,'bo',ms=6)
#            plt.plot(x_old,y_old,'wo',ms=4)
#            plt.plot(x_new,y_new,'bo')

            for i in xrange(len(x_old)):
                y_new[i] = y_new[i] - meanc
#            plt.plot(x_new,y_new,'ro')

            for i in xrange(len(x_old)):
                if y_new[i] < 0:
                   y_new[i] = abs(y_new[i])
                stacked_m.append(x_new[i])
                stacked_c.append(y_new[i])

#            plt.plot(x_new,y_new,'go')

        if plot : plt.plot(stacked_m,stacked_c,'go',label='Cluster galaxies')
        mags = []
        for mag in arange(min(limits[colour]['xlim']),max(limits[colour]['xlim'])+(dmag/2.0),dmag):

            d = {}
            d[mag] = []
            mags.append(d)
            
        for i in mags:
            lower = [gal > i.keys()[0] for gal in stacked_m]
            upper = [gal <= i.keys()[0]+dmag for gal in stacked_m]
            test = [lower[j]*upper[j] == 1 for j in xrange(len(lower))]
            indices = compress([val == True for val in test],range(len(test)))
            if len(indices) > 0:
                for gal in indices:
                    i[i.keys()[0]].append(stacked_c[gal])

        smooth_m = []
        smooth_c = []
        
        for i in mags:
            if len(i[i.keys()[0]]) > 0:
                so = sort(i[i.keys()[0]])
                if len(so) > 5:
                    so = so[-5:]
                    smooth_m.append(i.keys()[0])
                    smooth_c.append(median(so))
                else:
                    smooth_m.append(i.keys()[0])
                    smooth_c.append(mean(so))

        nanfilter = numpy.isnan(smooth_c)
        nanfilter = (nanfilter == False)
        nanfilter2 = numpy.isnan(smooth_m)
        nanfilter2 = (nanfilter2 == False)
        filter = (nanfilter & nanfilter2)
        
        smooth_c = list(smooth_c[filter])
        smooth_m = list(smooth_m[filter])

        fit = numpy.poly1d(numpy.polyfit(smooth_m,smooth_c, deg=poly))
        fit_x = []
        fit_y = []
        for mag in arange(min(limits[colour]['xlim']),max(limits[colour]['xlim'])+(dmag/20.0),dmag/10.0):
            fit_x.append(mag)
            y = 0.0
            order = range(len(fit)+1)
            order.reverse()
            for i in order:
                y = y + (mag**(int(i)))*fit[i]

            fit_y.append(y)

        if plot:
            plt.plot(smooth_m,smooth_c,'b-',label='Binned maxwidth estimate')
            plt.plot(smooth_m,smooth_c,'bo',)
            plt.plot(fit_x,fit_y,'r-',label='Width function')
            ymax = max(smooth_c)+0.5*max(smooth_c)
            ymin = 0
            plt.legend()
            magnitude = colour[-1]+' mag'
            plt.xlabel(magnitude)
            plt.ylabel('%s width'%(colour))
            plt.title('Polynomial order %d fitted width function in %s-detected clusters'%(poly,colour))
            plt.ylim(limits[colour]['ylim'])
            plt.ylim(0.0,0.6)
            plt.xlim(limits[colour]['xlim'])
            plt.show()
        

        return fit

    def print_sel(self):
        selection = self.selection
        colours = ['gr','ri','iz','zy']
        m = '_slope'
        c = '_intercept'
        w = '_width'
        pt = self.selection['detection']['probthresh']
        rc = self.selection['detection']['rho_crit']

        try:
            pt = self.selection['overrides']['probthresh']
        except:
            pass

        try:
            rc = self.selection['overrides']['rho_crit']
        except:
            pass
            
        mag = selection['defaults']['magnitude']
        string = ''
        for col in colours:
            string = string + '%s:m=%3.3f,c=%3.3f+/%3.1f\n'%(str(col),selection['colours'][col]['slope'],selection['colours'][col]['intercept'],selection['colours'][col]['width'])
        string = string + '%s<%4.1f\n'%(mag,selection['limits']['magnitudes'][mag]['limit'])
        string = string + 'z<%4.2f\n'%(selection['limits']['magnitudes']['zmax'])
        string = string + 'rho_crit=%4.1f\n'%(rc)
        string = string + 'probthresh=%4.4f'%(pt)
        return string



    def width_evol(self,plot=True,colour='default',perfect=False,passback=False,style='k.',compare=False,metric='quality',thresh=0.5,raw=False,gridsize=50,fit=True,contour=False,nbins=20,transpose=False,bootstrap=False,poly=2):


        if colour == 'default':
            try:
                colour = self.selection['defaults']['colour']
            except:
                colour = 'gr'
        cols = []
        widths = []
        for cluster in self.clusters:
            grad,cl_col,cl_width = cluster.cmr(fit=True,colour=colour)
            cols.append(cl_col)
            widths.append(cl_width)
            
        if plot:
            plt.plot(cols,widths,style)
            plt.ylabel('%s colour @ %s=20'%(colour,colour[-1]))
            plt.ylabel('width in %s'%(colour))
            plt.xlim(self.selection['limits']['cmr_range'][colour][0],self.selection['limits']['cmr_range'][colour][1])
            


    def slope_evol(self,plot=True,colour='default',perfect=False,passback=False,style='k.',compare=False,metric='quality',thresh=0.5,raw=False,gridsize=50,fit=True,contour=False,nbins=20,transpose=False,bootstrap=False,poly=2):
        """\
        Plots the colour normalisation-slope relation

        """
        if colour == 'default':
            try:
                colour = self.selection['defaults']['colour']
            except:
                colour = 'gr'


        import numpy
        if passback == False and plot == True:
            fig = plt.figure()
            sp = plt.subplot(111,aspect='equal')

        limits = self.limits
        magnitude = colour[-1]+'mag'
        col = 'galaxy.'+colour
        mag = 'galaxy.'+magnitude
        clusters = self.clusters
        if perfect:
            clusters = self.perfect_clusters
        col_m20 = []
        slope = []
        cm20_err = []
        slope_err = []
        for cluster in clusters:
            if bootstrap:
                grad,wth,cm20,gerr,werr,cerr = cluster.cmr(fit=True,bootstrap=True,error=True,colour=colour)
                slope_err.append(gerr)
                cm20_err.append(cerr)
                    
            else:
                grad,wth,cm20 = cluster.cmr(fit=True,colour=colour)
            col_m20.append(cm20)
            slope.append(grad)


        _fit = numpy.poly1d(numpy.polyfit(slope,col_m20,deg=1))
        
        if transpose:
            fit2 = numpy.poly1d(numpy.polyfit(col_m20, slope, deg=poly))
        else:
            fit2 = numpy.poly1d(numpy.polyfit(slope,col_m20, deg=poly))

        if raw:
            return [_fit,col_m20]
        fit_x,fit_y = [],[]
        fit_x.append(min(slope))
        fit_x.append(max(slope))
        fit_y.append((_fit[1]*fit_x[0])+_fit[0])
        fit_y.append((_fit[1]*fit_x[1])+_fit[0])
        
        if plot:
            plt.subplots_adjust(left=0.09,right=0.90,top=0.94,bottom=0.08,hspace=0.01)
#            plt.plot(slope,col_m20,style)
#            plt.plot(fit_x,fit_y,style[0]+':')
            if 1:
                x,y = slope,col_m20
                if contour:
                    if transpose:
                        plt.hexbin(y,x,gridsize=gridsize,bins='log',cmap=plt.cm.jet,antialiased=True,label='#clusters')
                    else:
                        plt.hexbin(x,y,gridsize=gridsize,bins='log',cmap=plt.cm.jet,antialiased=True,label='#clusters')
                    cb = plt.colorbar()
                    cb.set_label('log10(#clusters)')
                else:
                    if transpose:
                        if bootstrap:
                            plt.errorbar(y,x,xerr=cm20_err,yerr=slope_err,color='k',marker='o',ls='None')
                        else:
                            plt.plot(y,x,'ko')
                    else:
                        if bootstrap:
                            plt.errorbar(y,x,yerr=cm20_err,xerr=slope_err,color='k',marker='o',ls='None')
                        else:
                            plt.plot(x,y,'ko')

                if fit:
                    junk = """\
                    _x = []
                    _y = []
                    _xerr = []
                    _xwt = []
                    import numpy
                    dy = (max(y)-min(y))/float(nbins)
                    bins = arange(min(y),max(y),dy)
                    x = numpy.array(x)
                    y = numpy.array(y)
                    for bin in bins:
                        constraint = (y >= bin) & (y < bin+dy)
                        _selx = x[constraint]
                        _sely = y[constraint]
                        if len(_selx) > 0:
                            _y.append(bin+dy/2.0)
                            _x.append(median(_selx))
                            _xerr.append(std(_x))
                            _xwt.append(float(len(_x))/float(len(x)))
                        
                    _fit2 = numpy.poly1d(numpy.polyfit(_x,_y,deg=1))

                    import scipy
                    fitfunc = lambda p, x: (p[0]*x)+p[1]
                    errfunc = lambda p, x, y: fitfunc(p,x)-y
                    _x = numpy.array(_x)
                    _y = numpy.array(_y)
                    _thefit = scipy.optimize.leastsq(errfunc,(-60.0,-2.0),args=(_x,_y),full_output=0)
#                    print _thefit
                    _dx = max(x)-min(x)
                    _xrange = arange(min(x)-2.0*_dx,max(x)+2.0*_dx,0.2)
                    if transpose:
                        plt.plot(_xrange,fitfunc(_thefit[0],_xrange),'r-')
                    else:
                        plt.plot(fitfunc(_thefit[0],_xrange),_xrange,'r-')
                    _thefit = numpy.poly1d(_thefit[0])
                    """
                    
                    if transpose:
                        xmin,xmax = plt.xlim()
                        cm20range = arange(xmin,xmax,abs(xmax-xmin)/50.0)
#                        plt.plot(cm20range,fit2(cm20range),'r-')
                        plt.plot(cm20range,fit2(cm20range),'r-')
                        
                    
                    
                    style = 'w-'
                    if not contour:
                        style = 'k-'
#                    plt.plot([min(x),max(x)],[_fit2(min(x)),_fit2(max(x))],style,label='All clusters')
#                    plt.errorbar(_x,_y,xerr=_xerr,marker='o',color='r')
                if compare:
                    fit2 = self.slope_evol(plot=False,colour=colour,perfect=True,passback=False)[0]
                    plt.plot([min(x),max(x)],[fit2(min(x)),fit2(max(x))],'w--',label='Perfect,y=%s'%(str(fit2)[2:]))
                    import science
                    sci = science.science()
                    mets = sci.metrics(self)
                    filtered = sci.filter(self,metric.lower(),thresh,1.1)
                    fit3 = filtered.slope_evol(plot=False,colour=colour,perfect=False,passback=False)[0]
                    plt.plot([min(x),max(x)],[fit3(min(x)),fit3(max(x))],'w:',label='%s>%s,y=%s'%(metric[0].upper(),str(thresh),str(fit3)[2:]))


                try:
                    if transpose:
                        plt.xlim(self.selection['limits']['cmr_range'][colour][0],self.selection['limits']['cmr_range'][colour][1])
                    else:
                        plt.ylim(self.selection['limits']['cmr_range'][colour][0],self.selection['limits']['cmr_range'][colour][1])

                except:
                    if transpose:
                        plt.xlim(min(y),max(y))
                    else:
                        plt.ylim(min(y),max(y))

                if transpose:
                    plt.ylim(min(x),max(x))
                    ymin,ymax = plt.xlim()
                else:
                    plt.xlim(min(x),max(x))
                    ymin,ymax = plt.ylim()

                if transpose:
                    plt.plot([ymin,ymax],[0.0,0.0],'k--')
                    plt.ylabel('fitted slope')
                    plt.xlabel('%s colour @ %s=20'%(colour,magnitude))
                    plt.xlim(ymin,ymax)
                    plt.subplots_adjust(left=0.03,right=0.99)
                    
                else:
                    plt.plot([0.0,0.0],[ymin,ymax],'k--')                    
                    plt.xlabel('fitted slope')
                    plt.ylabel('%s colour @ %s=20'%(colour,magnitude))

                try:
                    plt.title('Evolution of slope in %s CMR, y=%s'%(colour,str(_fit)[2:]))
                except:
                    plt.title('Evolution of slope in %s CMR'%(colour))
                if 0:
                    lg = plt.legend(prop = matplotlib.font_manager.FontProperties(size='small'))
                    lg.get_frame().set_fill(False)
                    for t in lg.get_texts():
                        t.set_color('w')
                    frame = lg.get_frame()
                    frame.set_edgecolor('w')
                


        if fit:
            return [_fit,slope,col_m20,fit2]


        return [_fit,slope,col_m20]
    

    def wrap(self):
        from math import pi
        selection = self.selection
        try:
            for cl in self.clusters:
                cl.wrap()
        except:
            pass
        try:
            self.cuts['spatial']['x'] = self.cuts['spatial']['x'] - 2*pi
        except:
            pass
        try:
            for g in self.field_galaxies:
                g.wrap()
        except:
            pass
        try:
            for key in self.external_data.keys():
                self.external_data[key]['x'] = numpy.array(self.external_data[key]['x']) - 2*pi
        except:
            pass

        selection['limits']['xmin'] = selection['limits']['xmin'] - 2*pi
        selection['limits']['xmax'] = selection['limits']['xmax'] - 2*pi


        return

    def unwrap(self):
        from math import pi
        selection = self.selection
        try:
            for cl in self.clusters:
                cl.unwrap()
        except:
            pass

        try:
            self.cuts['spatial']['x'] = self.cuts['spatial']['x'] + 2*pi
        except:
            pass
        try:
            for g in self.field_galaxies:
                g.unwrap()
        except:
            pass
        try:
            for key in self.external_data.keys():
                self.external_data[key]['x'] = numpy.array(self.external_data[key]['x']) + 2*pi
        except:
            pass
        
        selection['limits']['xmin'] = selection['limits']['xmin'] + 2*pi
        selection['limits']['xmax'] = selection['limits']['xmax'] + 2*pi

        return

    def plot_selection(self,colour='default',sources=False,excluded=False,data=None,passback=False,style=['-',':'],plot_colour='r',verbose=False):
        #for a selection function with a slope, the intercept is defined as col when mag=20
        if colour == 'default':
            try:
                colour = self.selection['defaults']['colour']
            except:
                colour = 'gr'

        from clusters import width
        if not passback:
            plt.figure()
        selection = self.selection
        limits = self.limits

        if not self.clusters:
            self.clusters = []

        if selection['limits'].has_key('cmr_range') and len(self.clusters) > 0:
            lims = selection['limits']['cmr_range']
            try:
                mag = self.selection['colours'][colour]['magnitude']
            except:
                mag = colour[-1]
            galmags = []
            for c in self.clusters:
                galmags = galmags + [g.magnitudes[mag] for g in c.galaxies]
            xlim = [min(galmags)-1.0,max(max(galmags)+1.0,24.0)]
            limits = {colour:{'ylim':[lims[colour][0],lims[colour][1]],'xlim':xlim}}
        

        if selection['colours'][colour]['fitfuncs'].has_key('gradient_function'):
            fit = selection['colours'][colour]['fitfuncs']['gradient_function']
            gradient = (selection['colours'][colour]['intercept']-fit[0])/fit[1]
            print '\t#detections.plot_selection(): gradient=%f from intercept=%f'%(gradient,selection['colours'][colour]['intercept'])
        else:
            gradient = selection['colours'][colour]['slope']
            print '\t#detections.plot_selection(): Using simple gradient %f'%(gradient)

        intercept = selection['colours'][colour]['intercept']
        if gradient != 0:
            intercept = selection['colours'][colour]['intercept']-(gradient*20.0)

        mag = colour[-1]+'mag'
        if excluded:
            cols = []
            mags = []
            if not data:
                from dataIO import getData
                obs = False
                if self.selection.has_key('observed_data'):
                    obs = True
                try:
                    data = self.source_data
                    ll = len(data['spatial']['x'])
                except:
                    data = getData(self.selection['io']['sourcefile'],obs=obs,dataset=self.selection['io']['dataset'])
                try:
                    self.source_data = data
                except:
                    pass

            constraint = (data['rmag'] < selection['rlimit']) & (data[colour] >= limits[colour]['ylim'][0]) & (data[colour] <= limits[colour]['ylim'][1])
            mags = data[mag][constraint]
            cols = data[colour][constraint]
            plt.hexbin(mags,cols,gridsize=350,bins='log',cmap=plt.cm.binary,antialiased=True)
            if selection.has_key('cmr_range') and len(self.clusters) > 0:
                xlim = [min(mags)-1.0,max(max(mags)+1.0,24.0)]
                limits = {colour:{'ylim':[lims[colour][0],lims[colour][1]],'xlim':xlim}}    


        fit_x,fit_y = [],[]
        fit_x.append(min(limits[colour]['xlim']))
        fit_x.append(max(limits[colour]['xlim']))
        fit_y.append((gradient*fit_x[0])+intercept)
        fit_y.append((gradient*fit_x[1])+intercept)
        if not sources:
            plt.plot(fit_x,fit_y,color=plot_colour,ls=style[0])
        else:
            plt.plot(fit_x,fit_y,'k-')

        dmag = 0.1
        fit_x,fit_y,fit_y2 = [],[],[]
        for magnit in arange(min(limits[colour]['xlim']),max(limits[colour]['xlim'])+(dmag/20.0),dmag/10.0):
                fit_x.append(magnit)
                fit_y.append((gradient*fit_x[-1])+intercept+width(selection,colour,fit_x[-1],{colour:0.0}))
                fit_y2.append((gradient*fit_x[-1])+intercept-width(selection,colour,fit_x[-1],{colour:0.0}))
        if not sources:
            plt.plot(fit_x,fit_y,color=plot_colour,ls=style[1])
            plt.plot(fit_x,fit_y2,color=plot_colour,ls=style[1])
        else:
            plt.plot(fit_x,fit_y,'k:')
            plt.plot(fit_x,fit_y2,'k:')
        

        mags = []
        cols = []
        if sources:
            try:
                for i in xrange(len(self.source_galaxies['x'])):
                    mags.append(self.source_galaxies[mag][i])
                    cols.append(self.source_galaxies[colour][i])
            except:
                pass
        plt.plot(mags,cols,'b.')
        plt.xlim(limits[colour]['xlim'])
        plt.ylim(limits[colour]['ylim'])
        plt.xlabel(mag)
        plt.ylabel(colour)
#        plt.show()



    def cmr(self,colour='default',plot=False,fit=False,plotfits=False,bcgfit=False,override_limits=False,perfect=False,widthfit=False,nbins=40,object='clusters',contour=False):
        from numpy import polyfit
        if colour == 'default':
            colour = self.selection['defaults']['colour']
            magnitude = self.selection['colours'][colour]['magnitude']

        else:
            magnitude = self.selection['colours'][colour]['magnitude']
            
        if plot or plotfits:
            plt.figure()
        global pi
        import numpy,scipy
        limits = self.limits

        if 'gal' in object:
            print '\tChecking CMR of galaxies:'
        

        if override_limits: limits = {}
        limits['default'] = {'xlim':[14,25],'ylim':[-0.2,2.0]}
        clusters = self.clusters
        if clusters == None: return [0.,0.]
        if clusters == []: return [0.,0.]
        if perfect: clusters = self.perfect_clusters
        if plotfits:
            plot = True
            fit = True
        allm = []
        allc = []
        bcgm = []
        bcgc = []
        allfits = []
        for cluster in clusters:
            m = []
            c = []
#            col = 'galaxy.'+colour
#            mag = 'galaxy.'+magnitude
            for galaxy in cluster.galaxies:
                m.append(galaxy.magnitudes[magnitude])
                c.append(galaxy.colours[colour])
#            col = 'cluster.bcg.'+colour
#            mag = 'cluster.bcg.'+magnitude
            bcgm.append(cluster.bcg.magnitudes[magnitude])
            bcgc.append(cluster.bcg.colours[colour])
            allm = allm + m
            allc = allc + c
            
            m = numpy.array(m)
            c = numpy.array(c)

            nanfilter = numpy.isnan(c)
            nanfilter = (nanfilter == False)
            nanfilter2 = numpy.isnan(m)
            nanfilter2 = (nanfilter2 == False)
            filter = (nanfilter & nanfilter2)

            c = list(c[filter])
            m = list(m[filter])

            if len(m) >= 2:
                grad,int = polyfit(m,c,1)
                allfits.append(grad)


        if widthfit:
            import math
# need to find out if there was a gradient set to the previous colour...            
            rotate = False
            colour_list = self.selection['detection']['sequence']
            if colour == colour_list[0]:
                prev_colour = colour
                if self.selection['colours'][colour]['slope'] != 0: 
                    rotate = True
                    print '\tUsing this (%s) colour gradient to rotate data'%(colour)
            else:
                prev_colour = ''
                for i in xrange(len(colour_list)):
                    if colour == colour_list[i]:
                        try:
                            prev_colour = colour_list[i-1]
                            if self.selection['colours'][prev_colour]['slope'] != 0: 
                                rotate = True
                                print '\tUsing previous (%s) colour gradient to rotate data'%(prev_colour)
                        except:
                            print '\tNot expecting anything to happen here'
                        #first colour, so ignore...

            col = []
            mag = []
            bcgcol = []
            bcgmag = []
            parent = 'cluster.'
            source = 'bcg.'

            c = 'cluster.bcg.'+colour
            m = 'cluster.bcg.'+magnitude
            for cluster in clusters:
                    bcgcol.append(cluster.bcg.colours[colour])
                    bcgmag.append(cluster.bcg.magnitudes[magnitude])

            c = 'galaxy.'+colour
            m = 'galaxy.'+magnitude
            for cluster in clusters:
                _col,_mag = cluster.get_colours_mags(colour,magnitude)
                col = list(col) + list(_col)
                mag = list(mag) + list(_mag)


            mag = numpy.array(mag)
            col = numpy.array(col)
            if 1:
#                print len(mag),mag,len(col),col
                nanfilter1 = numpy.isnan(col)
                nanfilter1 = (nanfilter1 == False)
#                print nanfilter1
                nanfilter2 = numpy.isnan(mag)
                nanfilter2 = (nanfilter2 == False)
#                print nanfilter2
#                print len(nanfilter1), len(nanfilter2)
                nanfilter = (nanfilter1 & nanfilter2)

                col = list(col[nanfilter])
                mag = list(mag[nanfilter])
                
            if len(mag) >= 2:
#                print type(m),m,type(c),c
                grad,int = polyfit(mag,col,1)
                allfits.append(grad)


            if False in filter:
                allm = mag
                allc = col

            if len(col) == 0:
                return 0.0,0.0

            if rotate and self.selection['colours'][prev_colour]['fitfuncs'].has_key('gradient_function'):
                print '\t#cmr: detected %s gradient function - fitting slope from current clusters'%(prev_colour)
                self.selection['colours'][prev_colour]['slope'] = average(self.cmr(fit=True,colour=prev_colour))

            if rotate:
                gradient = self.selection['colours'][prev_colour]['slope']
                meanmag = mean(allm)
                meanc = mean(allc)
                y0 = 0.0
                x0 = (-1.0*meanc/gradient)+meanmag
                r = []
                angle = []
                x_old = []
                y_old = []
                bcgx_old = []
                bcgy_old = []
                bcg_r = []
                bcg_ang = []
                for i in xrange(len(allm)):
                    allm[i] = allm[i] - x0
                    allc[i] = allc[i] - y0
                    r.append(((allm[i]**2)+(allc[i]**2))**0.5)
                    x_old.append(allm[i]/r[i])
                    y_old.append(allc[i]/r[i])
                    angle.append(math.atan2(allc[i],allm[i]))
                for i in xrange(len(bcgm)):
                    bcgm[i] = bcgm[i] - x0
                    bcgc[i] = bcgc[i] - y0
                    bcg_r.append(((bcgm[i]**2)+(bcgc[i]**2))**0.5)
                    bcgx_old.append(bcgm[i]/bcg_r[i])
                    bcgy_old.append(bcgc[i]/bcg_r[i])
                    bcg_ang.append(math.atan2(bcgc[i],bcgm[i]))

                fit_angle = math.atan2(-1.0*gradient,-1.0)
                #need to rotate by 2pi - fit_angle
                d_theta = math.pi - fit_angle
                print '\tRotating by %f'%(d_theta)
                x_new = []
                y_new = []
                bcgx_new,bcgy_new = [],[]
                for i in xrange(len(x_old)):
                    angle[i] = angle[i] + d_theta
                    x_new.append(math.cos(angle[i]))
                    y_new.append(math.sin(angle[i]))

                for i in xrange(len(bcgx_old)):
                    bcg_ang[i] = bcg_ang[i] + d_theta
                    bcgx_new.append(math.cos(bcg_ang[i]))
                    bcgy_new.append(math.sin(bcg_ang[i]))

                for i in xrange(len(x_old)):
                    x_new[i] = r[i]*x_new[i]
                    y_new[i] = r[i]*y_new[i]
                    x_old[i] = r[i]*x_old[i]
                    y_old[i] = r[i]*y_old[i]

                    x_new[i] = x_new[i] + x0 
                    y_new[i] = y_new[i] + y0 + meanc
                    x_old[i] = x_old[i] + x0
                    y_old[i] = y_old[i] + y0
                
                for i in xrange(len(bcgx_old)):
                    bcgx_new[i] = bcg_r[i]*bcgx_new[i]
                    bcgy_new[i] = bcg_r[i]*bcgy_new[i]
                    bcgx_old[i] = bcg_r[i]*bcgx_old[i]
                    bcgy_old[i] = bcg_r[i]*bcgy_old[i]

                    bcgx_new[i] = bcgx_new[i] + x0
                    bcgy_new[i] = bcgy_new[i] + y0 + meanc
                    bcgx_old[i] = bcgx_old[i] + x0
                    bcgy_old[i] = bcgy_old[i] + y0

                
                if plot and widthfit == False:
                    plt.plot(x_old,y_old,'bo',ms=6)
                    plt.plot(x_old,y_old,'wo',ms=4)
                    plt.plot(bcgx_old,bcgy_old,'ro',ms=6)
                    plt.plot(bcgx_old,bcgy_old,'wo',ms=4)

                    plt.plot(x_new,y_new,'bo')
                    plt.plot(bcgx_new,bcgy_new,'ro')
                    x = arange(-10,30)
                    m_old = [x_old[i] for i in xrange(len(x_old))]
                    c_old = [y_old[i] for i in xrange(len(x_old))]
                    m_new = [x_new[i] for i in xrange(len(x_old))]
                    c_new = [y_new[i] for i in xrange(len(x_old))]
                    meanmag = mean(m_old)
                    meanc = mean(c_old)
                    m,p = polyfit(m_old,c_old,1)
                    y = (m*x)+meanc
                    x = x + meanmag
                    plt.plot(x,y,'-k')

                    x = arange(-10,30)
                    meanmag = mean(m_new)
                    meanc = mean(c_new)
                    m,p = polyfit(m_new,c_new,1)
                    y = (m*x)+meanc
                    x = x + meanmag
                    plt.plot(x,y,':k')

                    plt.ylim(limits[colour]['ylim'])
                    plt.xlim(limits[colour]['xlim'])
#                    plt.show()

                col,bcgcol = y_new,bcgy_new
                
            from scipy.optimize import leastsq
            import scipy
            import copy
            fitfunc = lambda p, x: p[0]*scipy.exp(-(x-p[1])**2/(2.0*p[2]**2))
            errfunc = lambda p, x, y: fitfunc(p,x)-y
            med = median(col)
            sdev = std(col)

            try:
                m1,m2 = self.limits[colour]['ylim'][1],self.limits[colour]['ylim'][0]
            except:
                this_mag = self.selection['colours'][colour]['magnitude']
                try:
                    m1,m2 = self.selection['limits']['magnitudes'][this_mag]['limit'],13.0
                except:
                    m1,m2 = 24.0,13.0
                
            delta_mag = m1-m2
            dbin = (delta_mag)/nbins
            bins = arange(m2,m1+dbin,dbin)

            outlier_metric = numpy.array([abs(gal_col/sdev) < 6.0 for gal_col in col])
            tmp_col = numpy.array(col)
            to_fit = list(tmp_col[outlier_metric])
            col_backup = deepcopy(col)
            if len(to_fit) != 0:
                col = to_fit

            #check, in case some colours are NaN (permit_dropouts=True)
            

            med = median(col)
            sdev = std(col)

#            limit_range = self.limits[colour]['ylim'][1]-self.limits[colour]['ylim'][0]
            limit_range = m1-m2

            if countOf(numpy.isnan(col),True) == len(col):
                x_gals_true = bins
                y_gals_true = numpy.array([0 for _tmp in bins])
            else:
                _col = col
                if True in numpy.isnan(col):
                    _col = []
                    for cc in col:
                        if not numpy.isnan(cc):
                            _col.append(cc)
                y_gals_true,x_gals_true = numpy.histogram(_col,bins=bins)
                if len(y_gals_true) != len(x_gals_true):
                    y_gals_true = list(y_gals_true)+[0]

                
            if countOf(numpy.isnan(bcgcol),True) == len(bcgcol):
                x_bcg_true = bins
                y_bcg_true = numpy.array([0 for _tmp in bins])
            else:
                _bcgcol = bcgcol
                if True in numpy.isnan(bcgcol):
                    _bcgcol = []
                    for cc in bcgcol:
                        if not numpy.isnan(cc):
                            _bcgcol.append(cc)
                y_bcg_true,x_bcg_true = numpy.histogram(_bcgcol,bins=bins)
                if len(y_bcg_true) != len(x_bcg_true):
                    y_bcg_true = list(y_bcg_true)+[0]


                 
#            try:
#               
#            except:
#                y_gals_true,x_gals_true = scipy.histogram(col,bins=bins)
#                y_bcg_true,x_bcg_true = scipy.histogram(bcgcol,bins=bins)

#            y_bcg_true = list(y_bcg_true)+[0]
            
                
            p0_gals = [max(y_gals_true),median(col),std(col)]
            p0_bcg = [max(y_bcg_true),median(bcgcol),std(bcgcol)]
#            print p0_gals
            fit = scipy.optimize.leastsq(errfunc,p0_gals,args=(x_gals_true,y_gals_true),full_output=1)
            pgals = fit[0]
            gal_status = fit[-1]
            message = fit[-2]

#            pbcg,bcg_status = scipy.optimize.leastsq(errfunc,p0_bcg,args=(x_bcg_true,y_bcg_true))
            fit_bcg = scipy.optimize.leastsq(errfunc,p0_bcg,args=(x_bcg_true,y_bcg_true))
            pbcg = fit_bcg[-2]
            bcg_status = fit_bcg[-2]
            bcg_message = fit_bcg[-2]
            passed = [1,2,3,4]
#            Disabled gaussian fit: cannot handle multiple gaussians in one window - fixing this would be good...
#            passed = [3.14159]

            p1 = copy.copy(pgals)
            kurtosis = p1[0]/p1[2]
            yatsdev = fitfunc(p1,p1[1]+p1[2])
            alt_kurtosis = (p1[0]-yatsdev)/p1[2]
            sd_coverage = 2.0*p1[2]/(m1-m2)
            col_coverage = (max(col)-min(col))/(m1-m2)
            simple_coverage = sdev/(m1-m2)

            if kurtosis < 50 and alt_kurtosis < 20 and col_coverage > 1:
                print '\tkurtosis,alt_kurtosis,sd_coverage,col_coverage,simple_coverage='
                print kurtosis,alt_kurtosis,sd_coverage,col_coverage,simple_coverage
                print '\tDetected poor gaussian form - overwriting sigma from %f to %f'%(sdev,sdev/20.0)
                message = 'maxfev'
                

            x_true,y_true = x_gals_true,y_gals_true
            if (not gal_status in passed) or (p1[2] < 0) or (p1[1] < m2) or (p1[1] > m1)\
                   or (str(p1[1]) == 'nan') or (p1[2] == 0) or ('maxfev' in message):
                p1 = copy.copy(pbcg)
#                x_true,y_true = x_bcg_true,y_bcg_true
                if True:
#                if not (bcg_status in passed) or p1[2] < 0:
                    p1[1] = med
                    p1[2] = sdev/20.0
                    print '\tSub-optimal fit to gaussian for galaxies and BCG galaxies'
                    print '\tSwitching to simple mean and stdev'
                    print '\tgalaxies : \t\t mean=%f sigma=%f'%(pgals[1],pgals[2])
                    print '\tBCGs     : \t\t mean=%f sigma=%f'%(pbcg[1],pbcg[2])
                    print '\tSimple   : \t\t mean=%f sigma=%f'%(med,sdev)
                    print '\tChosen   : \t\t mean=%f sigma=%f'%(p1[1],p1[2])
                    
            if plot:
                x_gals_model = scipy.linspace(min(x_gals_true),max(x_gals_true),50)
                y_gals_model = fitfunc(pgals,x_gals_model)
                x_bcg_model = scipy.linspace(min(x_bcg_true),max(x_bcg_true),50)
                y_bcg_model = fitfunc(pbcg,x_bcg_model)

                plt.plot(x_gals_true,y_gals_true,'k-',label='Cluster galaxies')
                plt.plot(x_gals_model,y_gals_model,'k:',label='Gaussian Fit')

                plt.plot(x_bcg_true,y_bcg_true,'b-',label='BCG Galaxies')
                plt.plot(x_bcg_model,y_bcg_model,'b:',label='Gaussian Fit')

                plt.ylabel('Frequency')
                plt.xlabel(colour)
                plt.legend(loc=1)
#                plt.show()

            return [p1[1],p1[2]]

        if bcgfit:
            plt.xlabel('r')
            plt.ylabel(colour)
            m,m2 = polyfit(bcgm,bcgc,1)
            x = arange(-100,100)
            meanmag = mean(bcgm)
            c = mean(bcgc)
            y = (m*x)+c
            x = x + meanmag
            plt.plot(bcgm,bcgc,'ro')
            plt.plot(x,y,':k')
            try:
                plt.ylim(limits[colour]['ylim'])
                plt.xlim(limits[colour]['xlim'])
            except:
                plt.ylim(limits['default']['ylim'])
                plt.xlim(limits['default']['xlim'])
                
#            plt.show()
            return [m,c]

        if plot:
            plt.plot(allm,allc,'bo')
            plt.plot(bcgm,bcgc,'ro')
            if not fit:
                plt.xlabel(magnitude)
                plt.ylabel(colour)
                plt.ylim(limits[colour]['ylim'])
                plt.xlim(limits[colour]['xlim'])
#                plt.show()

        if plotfits or (plot and fit):
            allfits = []
            plt.xlabel(magnitude)
            plt.ylabel(colour)
            x = arange(-10,30)
            col = 'galaxy.'+colour
            mag = 'galaxy.'+magnitude
            for cluster in clusters:
#                plt.ylim(limits[colour]['ylim'])
#                plt.xlim(limits[colour]['xlim'])
                x = arange(-10,30)
                m = [eval(mag) for galaxy in cluster.galaxies]
                c = [eval(col) for galaxy in cluster.galaxies]
                meanmag = mean(m)
                meanc = mean(c)
                m,p = polyfit(m,c,1)
                y = (m*x)+meanc
                x = x + meanmag
                allfits.append(m)
                plt.plot(x,y,':k')
            plt.ylim(limits[colour]['ylim'])
            plt.xlim(limits[colour]['xlim'])
#            plt.show()
            return allfits

        if fit:
            return allfits


    def info(self):
        print '\tInformation for this detection:'
        print '\t# galaxies in source selection: %d'%len(self.cuts['x'])
        try:
            print '\t# clusters: %d'%len(self.clusters)
            print '\t# galaxies in clusters: %d'%(len(self.cluster_galaxies))
        except:
            print '\t# clusters: %s'%(self.clusters)
            print '\t# galaxies in clusters: %s'%(self.cluster_galaxies)
        print '\tArea searched :%f sq deg'%(self.selection['detection']['detection_area'])
        print '\tSmallest cluster :%d'%(self.selection['detection']['smallest_cluster'])
        print '\tProbability threshold :%f'%(self.selection['detection']['probthresh'])
        print '\tCritical density :%f\n'%(self.selection['detection']['rho_crit'])

        print '\tColours:'
        for colour in self.selection['detection']['sequence']:
            print '\t%s: m=%f, c=%f, d(g-r)=%f'%(colour,self.selection['colours'][colour]['slope'],\
                                                 self.selection['colours'][colour]['intercept'],self.selection['colours'][colour]['width'])
        print '\tPositions:'
        print '\txmin=%f, xmax=%f'%(self.selection['limits']['xmin'],self.selection['limits']['xmax'])
        print '\tymin=%f, ymax=%f'%(self.selection['limits']['ymin'],self.selection['limits']['ymax'])
        try:
            print '\tzmin=%f, zmax=%f'%(self.selection['limits']['zmin'],self.selection['limits']['zmax'])
        except:
            pass


        try:
            print '\tPerformance:'
            print '\tCompleteness: %f'%(self.completeness())
            print '\tMean cluster purity: %f'%(mean(self.purity()))
            print '\tMean cluster quality: %f\n'%(mean(self.quality()))
        except:
            pass
        
    def performance(self,trainer,plot=False,label='',verbose=False,limits=None,bg='black',colours='default',units='default',grid=None,passback=False,size=6.0,r='auto',fill=False,observer_mode=True):
        """\
        detections.performance(self,trainer,plot=False) -> [f_det,f_field,cluster]
    
        \t Compare set of detected clusters to a training cluster

        trainer   : can be either a detections.detection object or a path
                    to a .cluster file
        cluster   : the cluster that matches the training cluster
        f_det     : fraction of galaxies detected that have been identified
                    in the training cluster (this fraction should be maximised):

                    #detected_cluster_galaxies/#training_cluster_members
                    
        f_field   : fraction of galaxies that were detected as cluster members
                    but were not indentified as members in the training cluster
                    (this fraction should be minimised):

                    #detected_galaxies_with_no_match_in_trainer/#detected_cluster_galaxies

                    caveat here is low completeness in the training cluster ie,
                    where it is assumed any galaxy not in the training cluster is
                    a field galaxy

        """
        import numpy
        if type(trainer) == type('a'):
            from dataIO import load
            try:
                train = load(trainer)
                if 'cluster' in trainer:
                    input_gals = train.galaxies
                    bcgname = train.bcg.name
                else:
                    input_gals = trainer.clusters[0].galaxies
                    bcgname = trainer.clusters[0].bcg.name
            except:
                print '\tFilename does not contain cluster information'

        else:
            input_gals = trainer.clusters[0].galaxies
            bcgname = trainer.clusters[0].bcg.name


        gx = [g.x for g in self.cluster_galaxies]
        if min(gx) < 0:
            for c in self.clusters:
                c.unwrap()


        gx = [g.x for g in trainer.cluster_galaxies]
        if min(gx) < 0:
            for c in trainer.clusters:
                c.unwrap()



        import numpy                
        gal_dna = numpy.array([hash('%1.8f %1.8f'%(float(g.x),float(g.y))) for g in self.cluster_galaxies]) 
        g = trainer.clusters[0].bcg
        trainer_bcg_dna = hash('%1.8f %1.8f'%(float(g.x),float(g.y)))
        t_dna = numpy.array([hash('%1.8f %1.8f'%(float(g.x),float(g.y))) for g in trainer.clusters[0].galaxies]) 

        if not trainer_bcg_dna in gal_dna:
            print '\tTraining cluster BCG not in cluster galaxies'
            from clusters import galaxies_within_r,clusters_within_r
            nearby = clusters_within_r(trainer.clusters[0],self.clusters,r=trainer.clusters[0].max_radius(),ranked=True)
            if len(nearby) == 0:
                nearby = galaxies_within_r(trainer.clusters[0],self.cluster_galaxies,r=trainer.clusters[0].max_radius())
                if len(nearby) > 0:
                    candidates = []
                    for n in nearby:
                        candidates.append(n.parent)
                    candidates = unique(candidates)
                else:
                    print 'No valid cluster found, returning nil points!'
                    return [0.0,1.0,trainer.clusters[0]]

            else:
                candidates = nearby
            fracs = []
            for c in candidates:
                c_dna = numpy.array([hash('%1.8f %1.8f'%(float(g.x),float(g.y))) for g in c.galaxies])
                members = [dna in t_dna for dna in c_dna]
                frac = float(countOf(members,True))/float(len(members))
                fracs.append(frac)
                candidate = candidates[indexOf(fracs,max(fracs))]
                
            cluster = candidate
            
        else:       
            cluster = self.cluster_galaxies[indexOf(gal_dna,trainer_bcg_dna)].parent
        detected_cluster_galaxies = []
        detected_field_galaxies = []
        detected_cluster_galaxies_trainer = []
        detected_field_galaxies_trainer = []


        cluster_dna = numpy.array([hash('%1.8f %1.8f'%(float(g.x),float(g.y))) for g in cluster.galaxies])
        trainer_dna = numpy.array([hash('%1.8f %1.8f'%(float(g.x),float(g.y))) for g in trainer.clusters[0].galaxies])

        reds = []
        greens = []

        for i in xrange(len(cluster_dna)):
            matches = numpy.where(trainer_dna == cluster_dna[i])[0]
            if len(matches) > 0:
                greens.append(cluster.galaxies[i])
            else:
                reds.append(cluster.galaxies[i])


        rdna = numpy.array([hash('%1.8f %1.8f'%(float(g.x),float(g.y))) for g in reds])
        gdna  = numpy.array([hash('%1.8f %1.8f'%(float(g.x),float(g.y))) for g in greens])

        f_detected = float(len(greens))/float(len(trainer.clusters[0].galaxies))
        f_field = float(len(reds))/float(len(cluster.galaxies))
        if verbose:
            print  f_detected,f_field
        if plot:
            from clusters import galaxies_within_r
            if passback == False:
                fig = plt.figure()
                sp = plt.subplot(111, axisbg=bg,aspect='equal')
            from clusters import cluster as cl
            try:
                if grid != None:
                    source_grid = grid
                else:
                    try:
                        source_grid = self.voronoi_grid
                    except:
                        pass
                g0 = source_grid[0]
#                print 'Shortcut'
            except:
                centroid = galaxy(0,cluster.x,cluster.y,dummy=True)
                dummy = cl([centroid],0,dummy=True)
                from tools import objectify
                cuts = deepcopy(self.cuts)
                if min(cuts['x']) < 0:
                    from math import pi
                    cuts['spatial']['x'] = cuts['spatial']['x'] + (2.0*pi)
            
                field = objectify(cuts)
                dummy.parent = self
                source_grid = dummy.local_density_field(field,r=cluster.max_radius()*5.0,plot=False,shape='square',wrap=False)
                self.voronoi_grid = source_grid

                
            sg = source_grid
            source_grid = galaxies_within_r(cluster,source_grid,r=cluster.max_radius()*5.0,shape='square')

            dna = numpy.array([hash('%1.8f %1.8f'%(float(g.x),float(g.y))) for g in source_grid])
            
            greens,reds,yellows,greys = [],[],[],[]
            
            for i in xrange(len(dna)):
                if dna[i] in gdna:
                    greens.append(source_grid[i])
                elif dna[i] in rdna:
                    reds.append(source_grid[i])
                elif dna[i] in trainer_dna:
                    yellows.append(source_grid[i])
                else:
                    greys.append(source_grid[i])

            from clusters import cluster as cl
            from detections import plot_galaxies

            cols = {'red':'r','green':'g','grey':'grey','yellow':'yellow'}
            
            fcol = {'red':None,'green':None,'grey':None,'yellow':None}
            if fill:
                fcol = cols
                if colours != 'default':
                    fcol = colours
                cols = {'red':'k','green':'k','grey':'k','yellow':'k'}
                
            elif colours != 'default':
                cols = colours
            colours = cols

            plot_galaxies(greys,colour=colours['grey'],passback=True,units=units,size=size,fill=fcol['grey'],observer_mode=observer_mode)
            plot_galaxies(reds,colour=colours['red'],passback=True,units=units,size=size,fill=fcol['red'],observer_mode=observer_mode)
            plot_galaxies(greens,colour=colours['green'],passback=True,units=units,size=size,fill=fcol['green'],observer_mode=observer_mode)

            if r == 'auto':
                rmax = trainer.clusters[0].max_radius()
                rmax = rmax * 2.0
            else:
                rmax = r
            xc,yc = trainer.clusters[0].x,trainer.clusters[0].y
            lims = {'limits':{'xmin':xc-rmax,'xmax':xc+rmax,'ymin':yc-rmax,'ymax':yc+rmax}}
            plot_galaxies(yellows,colour=colours['yellow'],passback=True,limits=lims,units=units,size=size,fill=fcol['yellow'],observer_mode=observer_mode)

            src = self
            if label == '':
                if src.selection['overrides'].has_key('rho_crit'):
                    rc = src.selection['overrides']['rho_crit']
                else:
                    rc = src.selection['detection']['rho_crit']

                if src.selection['overrides'].has_key('probthresh'):
                    pt = src.selection['overrides']['probthresh']
                else:
                    pt = src.selection['detection']['probthresh']
                label = 'P,rc=%1.3f,%3.1f f_det,f_field=%1.3f,%1.3f'%(pt,rc,f_detected,f_field)
            
            plt.title(label)
        return [f_detected,f_field,cluster]
#                    plt.show()
        

    def co_detect(self,colour,dc='auto',verbose=True,width=0.1,pp=False,ncpus=8,server=None,perfect=False,extras=None):
        print 'Passed to co_detect: server=%s'%(server)

        from copy import deepcopy,copy
        from clusters import cluster as new_cluster
        if pp:
            from parallel_detect import parallel_detect
            from parallel_vertices import deploy_vertices
            from parallel_qhull import qhull
            from dataIO import sel_edit
            from science import science
            sci = science()
            dependent_functions = (deploy_vertices,qhull,parallel_detect,sel_edit,sci.richness)
            modules = ('science','parallel_detect','clusters','numpy','parallel_qhull','parallel_vertices','vperc','time','detections')
            from dataIO import pp_run

        if self.clusters == None:
            print '\tCalculating firstpass clusters'
            clusters = parallel_detect(self.cuts,self.selection,plot=False)
            #don't want to add perfect clusters at this stage...
            self.addclusters(clusters,perfect=False)

        if dc == 'auto':
            dc = self.selection['limits']['cmr_range'][colour][-1]
            print 'Setting dc to %1.3f'%(dc)

        colours = self.selection['detection']['sequence']
        range = self.selection['limits']['cmr_range']

        m = self.selection['colours'][colour]['slope']
        try:
            test = len(m)
            m = (self.selection['colours'][colour]['intercept']-m[0])/m[1]
        except TypeError:
            m = self.selection['colours'][colour]['slope']
            
        core = deepcopy(self)
        
        c = self.selection['colours'][colour]['intercept']
        this_intercept,this_gradient = deepcopy(c),deepcopy(m)


        if verbose:
            print '\n\t####################################################################################'
            print '\t####################################################################################'
            print '\t                #co-detect in %s: c,m=%f,%f'%(colour,c,m)
            if self.selection['overrides'].has_key('parent_density'):
                print '\t                 [mean density from PARENT for snapshot]'
            print '\t####################################################################################'
            print '\t####################################################################################'
        cycle = iter(colours)
        c = ''
        lastcolour = False
        while c != colour:
            try:
                c = cycle.next()
            except StopIteration:
                if c != colour:
                    print '\tCannot find colour %s'%(colour)
                elif c == colours[-1]:
                    #this is the final band (eg iz), so cannot co-detect
                    return self,server
                else:
                    print '\tSomething unknown is wrong here...'
                    return self,server
        next_colour = cycle.next()
        mean,sig = self.cmr(widthfit=True,colour=next_colour)

        if self.selection['io']['dataset'] == 'snapshot':
            sig_set = 0.5
            print '#co-detect(): snapshot mode: overriding CMR %s colour search width from %f to %f'%(colour,sig,sig_set)
            sig = sig_set
            
        if verbose:
            print '\n\t####################################################################################'
            print '\t            #co-detect in %s: mean,sigma=%f,%f in next colour (%s)'%(colour,mean,sig,next_colour)

#        intercept_range = arange(mean-sig,mean+sig,dc)
        intercept_range = arange(max((mean-sig),core.selection['limits']['colours'][next_colour][0]),mean+sig,dc)      ##set blue limit of search to that of CMR exploration

        maxrange = arange(range[next_colour][0],range[next_colour][1],range[next_colour][2])
        if len(intercept_range) == 1: intercept_range = [mean]
        if len(intercept_range) == 0:
            print '\tNot expecting a zero-length intercept_range! setting to mean... (=%f)'%(mean)
            intercept_range = [mean]
        
        if len(intercept_range) > len(maxrange):
            if verbose:
                print '\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print '\t   #co-detect in %s: range too large from calculated spread - overriding to use'%(colour)
                print '\t                     default range for next colour (%s) with %d slices'%(next_colour,len(maxrange))
                print '\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            intercept_range = maxrange

        if verbose:

            print '\t            #co-detect in %s: creating %d slices'%(colour,len(intercept_range))
            print '\t####################################################################################'


#        print 'Current/next colour filter width=%f/%f'%(self.selection['colours'][colour]['width'],self.selection['colours'][next_colour]['width'])
#        raise SystemExit(0)

        if self.selection['colours'][next_colour]['width'] > 5.0:
            if self.selection['colours'][next_colour].has_key('trainer_width'):
                self.selection['colours'][next_colour]['width'] = self.selection['colours'][next_colour]['trainer_width']
                print 'Width of next colour filter has been set to its optimal from training cluster'
            else:
                self.selection['colours'][next_colour]['width'] = self.selection['colours'][colour]['width']
                print 'Width of next colour filter has been set that of the current colour'

        if self.selection.has_key('overrides'):
            if self.selection['overrides'].has_key(next_colour+'_width'):
                print '\tOverriding %s colour width from %f to %f'%(next_colour,self.selection['colours'][colour]['width'],self.selection['overrides'][next_colour+'_width'])
                self.selection['colours'][next_colour]['width'] = self.selection['overrides'][next_colour+'_width']
            elif self.selection['overrides'].has_key('fitted_width') and self.selection['colours'][next_colour].has_key('fitted_width'):
                print '\tOverriding %s colour width from %1.4f to %1.4f (fitted width)'%(next_colour,self.selection['colours'][colour]['width'],self.selection['colours'][next_colour]['fitted_width'])
                self.selection['colours'][next_colour]['width'] = self.selection['colours'][next_colour]['fitted_width']



        md_set = False
        if self.selection['overrides'].has_key('parent_density'):
            from detections import parent_mean_density
            print pp,server
            mean_densities,server = parent_mean_density(self,intercept_range,next_colour+'_intercept',job_server=server,ncpus=ncpus,pp=pp)
            md_set = True
            print '#detections.co_add(): recovered mean_densities (pp=%s)'%(pp)
            #print mean_densities

        if extras != None:
            if extras.has_key('form_detections'):
                from detections import detections as det
                base_selection = deepcopy(self.selection)
                dets = []
                for j in xrange(len(intercept_range)):
                    i = intercept_range[j]
                    base_selection['colours'][next_colour]['intercept'] = i
                    base_selection['colours'][next_colour]['intercept_index'] = j
                    w = base_selection['colours'][next_colour]['width']
                    if md_set:
                        base_selection['overrides']['mean_density'] = mean_densities[j]
                    dets.append(det(None,base_selection,dummy=True))
                parent_txt = ''
                if self.selection['overrides'].has_key('parent_density'):
                    parent_txt = ' (containing PARENT mean densities for snapshot mode)'
                print '#detections.co_detect(): Returning %d detections%s'%(len(dets),parent_txt)
                

                return dets
                

        source,server = datacuts.sel_edit2(self.selection,next_colour+'_intercept',intercept_range,pp=pp,perfect_clusters=False,ncpus=ncpus,server=server)
        if len(intercept_range) == 1: source = [source]

        for i in xrange(len(intercepts)):
            source[i].selection['colours'][next_colour]['intercept_index'] = i

        if self.selection['overrides'].has_key('parent_density'):
            for i in xrange(len(source)):
                source[i].selection['overrides']['mean_density'] = mean_densities[i]


            print '\t            #co-detect in %s: mean density override applied'%(colour)


        raise SystemExit(0)
        #### Why has this ^^ been put here, and why do we never get to this stage? I think accelerate() requests a return before we get down to here....


        co_detections = []
        if verbose:
            print '\n\t####################################################################################'
            print '\t            #co-detect in %s: now detecting clusters in %d slices'%(colour,len(intercept_range))



        if pp:
            run_function = 'parallel_detect'
            args = [{'source[i].cuts':None},{'source[i].selection':None}]
            co_detections,server = pp_run(source,run_function,args,dependent_functions,modules,verbose=False,ncpus=ncpus,server=server)
        else:
            from parallel_detect import parallel_detect
            for i in xrange(len(source)):
                co_detections.append(parallel_detect(source[i].cuts,source[i].selection))
        if verbose:
#            print '\tAdding these clusters to the slices'

            print '\t            #co-detect in c(%s)=%f: detected all clusters, adding to %d slices'%(colour,self.selection['colours'][colour]['intercept'],len(intercept_range))

        for i in xrange(len(intercept_range)):
            source[i].addclusters(co_detections[i],perfect=perfect)

        #now need to check that at least one cluster has been detected in at least one co-detection, otherwise we shouldn't call
        #cherrypick, and instead 'delete' the detection, by returning the detection class with no clusters (ie, this is a *conditional*
        #detection - a detection must be made in the original colour (eg g-r) AND the following colour (r-i).
        
        nclus = 0
        for det in source:
            nclus = max(nclus,len(det.clusters))

        if nclus == 0:
            best = copy(self)
            best.clusters = []
            if verbose:
                print '\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print '\t   #co-detect in %s: co-detection in %d slices returned zero clusters in all slices'%(colour,len(source))
                print '\t   producing and returning Null detection with no clusters (this is not the end of the world...)'
                print '\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        else:
            try:
                sci = science()
            except:
                try:
                    import science
                except ImportError:
                    print '\tCannot import science, not computing any more...'
                    return None,server
            try:
                sci = science.science()
            except:
                sci = science()
            if verbose:
                print '\t            #co-detect in %s: Cherrypicking clusters from %d slices'%(colour,len(intercept_range))
                print '\t####################################################################################'

            perfect_clusters = False
            best = sci.cherrypick(source,cm20_thresh=2.0*self.selection['colours'][colour]['width'],include_solo=True,colour=colour,perfect=perfect)
            if perfect:
                best,perfect_clusters = best
            else:
                best = best[0]

            
            if type(best) == type([1,2,3]) and len(best) > 0:
                from clusters import galaxy,cluster
                from dataIO import make_core
                new_det = detections(None,core.selection,dummy=True)
#                new_det = make_core(core.selection,perfect=False,dummy=True)
                new_det.addclusters(best,perfect=perfect_clusters)
                best = new_det

            elif type(best) == type([1,2,3]) and len(best) == 0:
                best = copy(core)
                best.clusters = []
                
                
                

                
        best.co_detections = True
        best.colour_colour = {}
        best.colour_colour['next'] = {'top':mean+sig,'middle':mean,'bottom':mean-sig,'name':next_colour}
        mean,sig = self.cmr(widthfit=True,colour=colour)
        best.colour_colour['this'] = {'top':mean+sig,'middle':mean,'bottom':mean-sig,'name':colour}
        best.cherry_flavour = 'inter'

        if verbose:
            print '\t\n'
            print '\t                #co-detect in %s: c,m=%f,%f   FINISHED'%(colour,this_intercept,this_gradient)
            print '\t####################################################################################'
            print '\t####################################################################################'
            print '\t\n'
            print '\t\n'
            print '\t\n'
            

        if pp:
            return best,server
        return best






def get_input(detection,perfect=False,data=None,pick_list=None):
    detection.get_input(data=data,pick_list=pick_list,perfect=perfect)
    return detection

def co_detect(detection,colour,dc=0.02,ncpus=8,width=0.1,pp_spawned=False,pp=False,perfect=False,verbose=True,pp_server=None,extras=None):
    print 'base co_detect: pp,server=%s,%s,%s'%(pp,pp_server,type(pp_server))
    cd = detection.co_detect(colour,dc=dc,verbose=verbose,width=width,pp=pp,ncpus=ncpus,server=pp_server,perfect=perfect,extras=extras)
    if pp:
        cd,server = cd[0],cd[1]              #why is this here?
        if pp_spawned == False: 
            #this process has not been spawned by pp_run: return the server (rather than "trash" it)
            return cd,server
    return cd

def plot_clusters(clusterlist,label='',colour=None,max_viz=False,limits=None,passback=False,centroids=False,size=4.0,fill=None,field_gals=None,naughty=False,interact=False,ext=True,rings=False,print_friendly=False,observer_mode=True,bcg_colour='yellow',lw=1.0):
    """\
        plot_clusters(clusterlist,spatialtest=False,label='',colour=None,max_viz=False,limits=None,passback=False,centroids=False,size=6.0,field_gals=None): -> void
    
        \t Plot an array of clusters.

        label     : Override label
        colour    : Override the plot colour [uses matplotlib standards - see plt.plot()]
        max_viz   : Override cluster colours and select from a colour table to increase contrast between clusters
        limits    : pass a selection class to define FOV limits
        passback  : Enable to tag this list onto an existing figure (be sure to specify
                    limits if this is the final set to plot
        bg        : override the background colour [default=black]
        centroids : Plot only the cluster centre positions
        size      : overide the default point size for the galaxies
        field_gals: see docstring for cluster.plot()
        interact  : query the data plot with the interaction class
        ext       : plot external data if available
        rings     : plot ring of radius r_max about cluster centroid
        lw        : override the VT linewidth
        """


# if spatial test, assume two clusters in clusterlist, need to shift one slightly and 
# use alt. plotting scheme


    from dataIO import Progress 

    if clusterlist == None:
        return

    if interact:
        try:
            import detections as dets
            container_class = dets.detections(None,None,dummy=True)
            container_class.selection = clusterlist[0].parent.selection
#            container_class = deepcopy(clusterlist[0].parent)
            container_class.clusters = clusterlist
            parents = [c.parent for c in clusterlist]
            gals = []
            for c in clusterlist:
                c.parent = container_class
                gals = gals + list(c.galaxies)
            container_class.cluster_galaxies = gals
#            container_class.addclusters(clusterlist)
            container_class.plot(interact=True,field=None,field_env=False,max_viz=max_viz,rings=rings,centroids=centroids,colour=colour,print_friendly=print_friendly,label=label,bcg_colour=bcg_colour,lw=lw)
            for c in clusterlist:
                c.parent = parents[indexOf(clusterlist,c)]
        except:
            print '\tProblem with producing container detections class - do these clusters belong to a detection?'
        return
    
    status = None
    if len(clusterlist) > 10:

        status = Progress(len(clusterlist),widget=False)


    if not passback:
        bg = 'black'
        if print_friendly:
            bg = 'white'
        fig = plt.figure()
        sp = plt.subplot(111,axisbg=bg,aspect=True)
    if max_viz:
        palatte = ['r','g','b','y','gray','w','magenta','cyan','orange','pink']
        paint_gun = iter(palatte)

    spatialtest = False
    if spatialtest and len(clusterlist) == 2:
        import copy
        newlist = []
        clusterlist[0].plot(passback=True,offset=clusterlist[0].x*0.001,centroid=centroids,galsize=size,fill=fill,ext=ext)
        clusterlist[1].plot(passback=True,centroid=centroids,galsize=size,fill=fill,ext=ext)


    else:
        for cluster in clusterlist:
            if max_viz:
                try:
                    colour = paint_gun.next()
                except:
                    paint_gun = iter(palatte)
                    colour = paint_gun.next()
            cluster.plot(passback=True,colour=colour,centroid=centroids,galsize=size,field_gals=field_gals,naughty=naughty,fill=fill,ext=ext,rings=rings,bcg_colour=bcg_colour)
            if status: status.update(0)

    if limits != None:
        lims = limits
        plt.xlim(lims['limits']['xmin'],lims['limits']['xmax'])
        plt.ylim(lims['limits']['ymin'],lims['limits']['ymax'])

    else:
        cx = [c.x for c in clusterlist]
        cy = [c.y for c in clusterlist]
        plt.xlim(min(cx),max(cx))
        plt.ylim(min(cy),max(cy))
    if not passback:
        plt.xlabel('RA')
        plt.ylabel('DEC')
    if label != '':
        plt.title(label)

    if observer_mode:
        x1,x2 = plt.xlim()
        plt.xlim(max(x1,x2),min(x1,x2))



#    plt.show()
        
def plot_galaxies(galaxylist,label=None,limits=None,vertices=True,passback=False,colour=None,style='-',bg='black',size=6.0,fill=None,verbose=False,ext_data=None,units='default',smart=False,ms='o',observer_mode=True,legend=None,lw=1.0):
    """\
        plot_galaxies(galaxylist,label='',limits=None,vertices=True,passback=False,colour=None,style='-',bg='black',size=6.0) -> void
    
        \t Plot an array of galaxies.

        label:    override label
        limits:   pass a selection class to define FOV limits
        vertices: Enable/disable plotting of galaxy vertices (False = faster plot)
        passback: Enable to tag this list onto an existing figure (be sure to specify
                  limits if this is the final set to plot
        colour:   override the plot colour [uses matplotlib standards - see plt.plot()]
        style:    override the vertex line style [uses matplotlib standards - see plt.plot()]
        bg:       override the background colour [default=black]
        size:     override the default point size for the galaxy
        legend:   set text for legend (will attempt to apply to a line, but else a point)
        lw:       define the linewidth of the vertices etc...
        """

    if verbose:
        from dataIO import Progress 
        status = Progress(len(galaxylist),widget=False)
        

    legend_text = [None for g in galaxylist]
    if len(legend_text) == 0:
        legend_text = [None]
    if legend != None:
        legend_text[0] = legend
    lt = legend_text
    legend_text = iter(legend_text)
        

    if not passback:
        fig = plt.figure()
        sp = plt.subplot(111, axisbg=bg,aspect='equal')

    if vertices:
        try:
            for galaxy in galaxylist:
#                galaxy.plot(passback=True,colour=colour,style=style,size=size,fill=fill,single_plot=False,ms=ms,legend=legend_text.next())
                galaxy.plot(passback=True,colour=colour,style=style,size=size,fill=fill,single_plot=False,ms=ms,legend=None,lw=lw)
                if verbose: status.update(0)
        except:
            #maybe no vertex info?
            pass
        g_x = [g.x for g in galaxylist]
        g_y = [g.y for g in galaxylist]
        if colour != None:
            col = colour
        else:
            col = 'b'

        if size > 0:
            plt.plot(g_x,g_y,color=col,marker=ms,ls='None',ms=size,label=legend_text.next())
            
    else:
        if not colour:
            col = 'r'
        else:
            col = colour
        g_x = []
        g_y = []
        for g in galaxylist:
            g_x.append(g.x)
            g_y.append(g.y)
        plt.plot(g_x,g_y,color=col,marker=ms,ls='None',ms=size,label=legend)

    if ext_data != None:
        ext_data.extdata_plot(passback=True,colour='red',lw=lw)

    if limits != None:
        lims = limits
        plt.xlim(lims['limits']['xmin'],lims['limits']['xmax'])
        plt.ylim(lims['limits']['ymin'],lims['limits']['ymax'])

    if units != 'default':
        from numpy import radians
        if units == 'degrees':
            xmin,xmax = plt.xlim()
            ymin,ymax = plt.ylim()
            dr = 0.1
            xvals = arange(degrees(xmin),degrees(xmax),dr)
            xrads = radians(xvals)
            if xvals[-1] < 0:
                from math import pi
                xvals = xvals + (360.0)
            newx = ['%2.1f'%(_x) for _x in xvals]

            init = newx[0].split('.')[0]+'.'
            test = [init in _val for _val in newx]
            if countOf(test,True) == len(test) and smart == True:            
                xv2 = [xv-xvals[0] for xv in xvals]
                xv2[0] = xvals[0]
                newx2 = ['+%2.1f'%(_x) for _x in xv2]
                newx2[0] = '%2.1f'%(xvals[0])
                newx = newx2

            yvals = arange(degrees(ymin),degrees(ymax),dr)
            yvals = yvals[1:]
            yrads = radians(yvals)
            newy = ['%2.1f'%(_y) for _y in yvals]


            plt.xticks(xrads,newx)
            plt.yticks(yrads,newy)

        elif units == 'None' or units == 'none':
#            print 'No units pls'
            frame1 = plt.gca()
#            frame1.axes.get_xaxis().set_visible(False)
#            frame1.axes.get_yaxis().set_visible(False)
            frame1.axes.get_xaxis().set_ticks([])
            frame1.axes.get_yaxis().set_ticks([])


            junk = """\
            for xlabel_i in frame1.axes.get_xticklabels():
                xlabel_i.set_visible(False)
                xlabel_i.set_fontsize(0.0)
            for xlabel_i in frame1.axes.get_yticklabels():
                xlabel_i.set_fontsize(0.0)
                xlabel_i.set_visible(False)
            for tick in frame1.axes.get_xticklines():
                tick.set_visible(False)
            for tick in frame1.axes.get_yticklines():
                tick.set_visible(False)
            """

    if observer_mode:
        x1,x2 = plt.xlim()
        plt.xlim(max(x1,x2),min(x1,x2))


    if not passback:
        plt.xlabel('RA')
        plt.ylabel('DEC')
        
        if label != None:
            plt.title(label)
    else:
        return


def common_cluster(detections,x,y,tolerance=1.0,plot=False,style='colour'):
    from clusters import galaxies_within_r
    from numpy import radians
    #assume tolerance in arcmins:
    tolerance = radians(tolerance/60.0)

    centroid = galaxy(-1,x,y)
    all_clusters = []
    commons = []
    for d in detections:
        d.compare('../bcg.dat')
        sel = d.selection
        this_list = []
        ints = {'gr':sel['gr_intercept'],'ri':sel['ri_intercept']}
        nearby = galaxies_within_r(centroid,d.cluster_galaxies,r=tolerance,shape='circle')

        if len(nearby) > 0:
            cluster_names = unique([g.parent.name for g in nearby])
            clusters = unique([g.parent for g in nearby])
            this_list = clusters
            all_clusters = all_clusters + list(clusters)

        commons.append([ints,this_list])


    
    if plot:
        import matplotlib
        extents = []
        for d in commons:
            clusters = d[1]
            this_det = []
            for c in clusters:
                this_det.append(c.max_radius())
            if len(this_det) > 0:
                extents.append(max(this_det))
        extents.sort()
        if len(extents) > 3:
            win = extents[-3]
        else:
            win = max(extents)
        sel = deepcopy(detections[0].selection)
        sel['limits']['xmin'] = x-win
        sel['limits']['xmax'] = x+win
        sel['limits']['ymin'] = y-win
        sel['limits']['ymax'] = y+win

        nplots = len(commons)

        lims = commons
        lims.reverse()
        iterator = iter(lims)
        d = iterator.next()
        while len(d[1]) == 0:
            d = iterator.next()

        lims.reverse()
        last_entry = indexOf(lims,d)
        nplots = last_entry + 1


        nrows,ncols = int(nplots**0.5)+1,1
        while nrows*ncols < nplots:
            ncols = ncols + 1
        
        if nplots < 5:
            nrows,ncols = nplots,1
        fig = plt.figure()
        col,row = 1,1
        i = 0
        for d in commons:
            if indexOf(commons,d) > nplots:
                break
            i = i + 1
            fig.add_subplot(nrows,ncols,i,axisbg='black',aspect='equal')
            if len(d[1]) > 0:
                plt.plot([x],[y],'ro')
                plot_clusters(d[1],limits=sel,passback=True)
                xlim = plt.gca().get_xlim()
                ylim = plt.gca().get_ylim()
                plt.plot([x],[y],'ro')
                plt.gca().set_xlim(xlim)
                plt.gca().set_ylim(ylim)

            plt.xticks([])
            plt.yticks([])

#            plt.title('ri=%1.3f,fd=%4.2f'%(d[0]['ri'],d[1][0].flux_density),size='small')
            if len(d[1]) > 0:
                if style == 'colour':
                    plt.title('gr=%1.3f,ri=%1.3f'%(d[0]['gr'],d[0]['ri']),size='small')   
                elif style == 'info':
                    plt.title('gr=%1.3f,ri=%1.3f,#=%d'%(d[0]['gr'],d[0]['ri'],len(d[1][0].galaxies)),size='small')   
                else:
                    plt.title('ri=%1.3f,fd=%4.2f,mfd=%1.4f'%(d[0]['ri'],d[1][0].flux_density/1E6,d[1][0].flux_density/(len(d[1][0].galaxies)*1E6)),size='small')
            else:
             plt.title('gr=%1.3f,ri=%1.3f'%(d[0]['gr'],d[0]['ri']),size='small')   
        plt.subplots_adjust(left=0.01,right=0.99,top=0.95,bottom=0.02,hspace=0.2,wspace=0.2)
    return commons
            
def darb(det,sampling_freq=10,save=False,area_only=False,healpix=False,nside=128):
#    print 'Starting sub_id %d'%(det.selection['detection']['subset_index'])
#    return [det,None]
    out = None

    if 1:
#    try:
        out = det.darb(sampling_freq=sampling_freq,status=False,save=save,area_only=area_only,healpix=healpix,nside=nside)

    else:
#    except:
        print 'Failed for sub_id %d'%(det.selection['detection']['subset_index'])
        return [None,None]

    if out == None:
        print '#darb(): I think this is an empty cell...'
        return [-1,None]

    hp = ''
    if healpix:
        hp = 'HEALPIX (nside=%s)'%(str(nside))
    print 'Finished sub_id %d, with %d samp freq, found %s area = %f'%(det.selection['detection']['subset_index'],sampling_freq,hp,det.darb_area)
    return [det,None]

def cluster_richness(det,bgcR=False,r80=False,bcg_centre=False,do_all=False):
    det.cluster_richness(pp=False,status_update=False,bgcR=bgcR,r80=r80,bcg_centre=bcg_centre,do_all=do_all)
    try:
        del det.subset_cuts
    except:
        pass
    return det

def cluster_simple_richness(det,bgcR=False,r80=False,bcg_centre=False,do_all=False,uber=False,rfrac=0.8):
#    print 'Uber mode =%s'%(str(uber))
    det.cluster_simple_richness(pp=False,status_update=False,bgcR=bgcR,r80=r80,bcg_centre=bcg_centre,do_all=do_all,uber=uber,rfrac=rfrac)
#    del det.src_data

#    try:
#        del det.subset_cuts
#    except:
#        pass
    return det




def retrieve_associates(det,keep_cuts=False):
    det.retrieve_associates(pp=False,status=False)
    if keep_cuts == False:
        try:
            del det.subset_cuts
        except:
            pass
    return det

def boundary_check(det,ra,dec,constraint=None,status=None):
    a = det.boundary_check(ra,dec,constraint=constraint,status=status)
    print 'Completed'
    return a


def parent_mean_density(source,variables,var_name,grid=None,job_server=None,ncpus=8,pp=False):
    """
    for the source (detection class) and a range of variables of name 'var_name', determine 
    the mean density from the parent catalogue and return these values for insertion as overrides

    example - you want to knowin advance the mean density of photometric cuts made scanning the CMR:
    variables = list [1.0,1.1,1.2,1.3...]
    var_name = 'gr_intercept'
    

    """
    #        0) duplicate selection, remove overrides and re-point the source/file, change xy limits [look for 'source_catalogue' override]
    #need to 1) generate array for cuts for full survey
    #        2) pass to detect_rapid with set_density=True
    #        3) add mean density overrides to each equivalent snapshot 
    import os
    from dataIO import setcfgProp, getcfgProp, getData, pp_run, sel_edit2
    from parallel_detect import parallel_detect
    from parallel_vertices import deploy_vertices
    from parallel_qhull import qhull
    from science import cherrypick
    from science import science
    import dataIO
    sci = science()
    if not source.selection['overrides'].has_key('source_catalogue'):
        print 'Requested mean densities from source catalogue, but have not'
        print """provided the source location. Please add 'source_catalogue' in .cfg file"""
        raise SystemExit(0)
    run_function = 'parallel_detect'
    dependent_functions = (deploy_vertices,qhull,parallel_detect,sci.richness,cherrypick,dataIO.Progress,sel_edit2)
    modules = ('science','parallel_detect','clusters','numpy','parallel_qhull','parallel_vertices','vperc','time','detections')
    args = [{'source[i].cuts':None},{'source[i].selection':None}]

    from copy import deepcopy
    dataset_src = source.selection['overrides']['source_catalogue']
    old_cfg_src = None
    if '/' in dataset_src:
        #check to see if it points to the same location as in orca.cfg, if not replace and swap back
        dataset_src_root = '/'.join(dataset_src.split('/')[:-1])+'/'
        dataset_src = dataset_src.split('/')[-1]
        dataset_src_root = os.path.abspath(dataset_src_root)
        if dataset_src_root[-1] != '/':
            dataset_src_root += '/'
        old_cfg_src = os.path.abspath(getcfgProp('orca','cfgpath',cfgloc='./'))
        if old_cfg_src[-1] != '/':
            old_cfg_src += '/'
        if old_cfg_src != dataset_src_root:
            #replace
            setcfgProp('orca','cfgpath',dataset_src_root,cfgloc='./')
    if '.cfg' in dataset_src:
        dataset_src = dataset_src.split('.cfg')[0]

    env_args = {}
    env_args['dc'] = source.selection['limits']['colours'].values()[0][-1]
    env_args['probthresh'] = source.selection['detection']['probthresh']
    env_args['rho_crit'] = source.selection['detection']['rho_crit']
    env_args['smallest_cluster'] = source.selection['detection']['smallest_cluster']
    env_args['halomass'] = source.selection['mock_data']['halomass']
    env_args['outlier_metric'] = source.selection['mock_data']['outlier_metric']
    env_args['metric'] = source.selection['detection']['metric']
    env_args['outlier_no_clobber'] = source.selection['mock_data']['no_clobber']
    
    _data,_junk,_junk,_range,_core = getData('',args=env_args,dataset=dataset_src,wrap=False,selection=None,subset_index=None,plot=False,grid=grid)
    sel_copy = _core.selection
    sel_copy['colours'] = deepcopy(source.selection['colours'])
    sel_copy['overrides'] = {}
    sels_parent,job_server = sel_edit2(sel_copy,var_name,numpy.array(variables),perfect_clusters=False,pp=pp,ncpus=ncpus,sourcedata=None,server=job_server,verbose=True)
    if type(sels_parent) != type([]):
        sels_parent = [sels_parent]

    sd_args = deepcopy(args)
    sd_args.append({'set_density':True})
    sd_args.append({'perfect':False})
    sd_args.append({'plot':False})
    sd_args.append({'verbose':False})
    
    mean_densities = []
    if pp:
        mean_densities,job_server = pp_run(sels_parent,run_function,sd_args,dependent_functions,modules,ncpus=ncpus,idle_delay=360.0,verbose=True,server=job_server)
    else:
        for s in sels_parent:
            mean_densities.append(parallel_detect(s.cuts,s.selection,plot=False,set_density=True,verbose=False))

    if old_cfg_src != None:
        if os.path.abspath('./')+'/' == old_cfg_src:
            old_cfg_src = './'
        setcfgProp('orca','cfgpath',old_cfg_src,cfgloc='./')

    print '#detections.parent_mean_density(): parent survey mean densities calculated'
    return mean_densities,job_server







def accelerate(detections,dc='auto',job_server=None,perfect=False,colour='auto',full_data=None,richness=True,ncpus=16):
#    from parallel_vertices import deploy_vertices
#    from parallel_qhull import qhull
    from parallel_detect import parallel_detect
    from science import cherrypick,science
    import dataIO
    from dataIO import save,load,pp_run,selection_fn,slice_data,make_mask,make_detections,getData,sel_edit,sel_edit2,Progress
    import fits_table as fits
    from detections import detections as det
    from numpy import unique
    run_function = 'co_detect'
    verbose = True
    extras = {'form_detections':True}

    sci = science()

    uber_det = detections[0].clone()
    uber_det.src_data = {}
    for c in uber_det.selection['detection']['sequence']:
        uber_det.selection['colours'][c]['width'] = 100.0
    if full_data:
        uber_det.src_data['total'] = full_data
    else:
        uber_det.get_input()
        uber_det.src_data['total'] = uber_dets.cuts
    

    try:
        permit_dropouts = detections[0].selection['io']['permit_dropouts']
    except KeyError:
        permit_dropouts = False

    if colour =='auto':
        colour = detections[0].selection['defaults']['colour']

    dataset = detections[0].selection['io']['dataset']

    if dc == 'auto':
        dc = detections[0].selection['limits']['cmr_range'][colour][-1]


    if permit_dropouts:
        first_stage = detections
        cl_members = {}
        cl_gals = {}
        cl_gals_obj = {}
        galID_lookup = {}
        cl_colourID_lookup = {}
        for det in detections:
            for cl in det.clusters:
                cl_members[cl.dna] = []
                if not cl.colour_dna() in cl_colourID_lookup.keys():
                    cl_colourID_lookup[cl.colour_dna()] = []
                cl_colourID_lookup[cl.colour_dna()].append([cl,det])

            if len(det.clusters) > 0:
                for g in det.cluster_galaxies:
                    cl_gals[g.extras['id']] = []
                    cl_gals_obj[g.extras['id']] = []
                
        for det in detections:
            for cl in det.clusters:
                cl_members[cl.dna] = cl_members[cl.dna] + list([g.extras['id'] for g in cl.galaxies])
            if len(det.clusters) > 0:
                for g in det.cluster_galaxies:
                    cl_gals[g.extras['id']] = cl_gals[g.extras['id']] + [g.parent.dna]
                    cl_gals_obj[g.extras['id']] = cl_gals_obj[g.extras['id']] + [g.parent]
                    galID_lookup[g.extras['id']] = g

                
        #sanity check (only if magmod is not being used, otherwise this gets triggered...):
        try:
            depth_mod = det.selection['overrides']['depth_mod'] > 0
        except:
            depth_mod = False            
        from numpy import unique
        if depth_mod == False:
            for key in cl_gals_obj:
                if len(unique([_cl.parent for _cl in cl_gals_obj[key]])) != len([_cl.parent for _cl in cl_gals_obj[key]]):
                    print 'Error - galaxyID %d has been assigned to >1 cluster in the same detection'%(key)
                    raise SystemExit(0)

        #    out = co_detect(detections[0],colour=colour,dc=dc,ncpus=ncpus,width=99,pp_spawned=False,pp=False,perfect=False,verbose=True,server=False,extras=extras)
    dependent_functions = (deploy_vertices,qhull,parallel_detect,sel_edit,sci.richness,cherrypick,dataIO.Progress,sel_edit2)
    modules = ('science','parallel_detect','clusters','numpy','parallel_qhull','parallel_vertices','vperc','time','detections')
    run_function = 'co_detect'
    pp_args = [{'source[i]':None},{'colour':colour},{'dc':dc},{'ncpus':ncpus},{'width':99},{'pp_spawned':True},{'pp':False},{'perfect':('mock' in dataset and perfect)},{'verbose':True},{'pp_server':False},{'extras':extras}]

    if 0:
        #to run with pp - experimental...
        pp_args = [{'source[i]':None},{'colour':colour},{'dc':dc},{'ncpus':ncpus},{'width':99},{'pp_spawned':True},{'pp':True},{'perfect':('mock' in dataset and perfect)},{'verbose':True},{'pp_server':'None'},{'extras':extras}]
                
    #pp_args = args

    if 0:
        _test = detections[21].clone(full=True)
#        detections[21].co_detect(colour,dc=dc,width=99,pp=False,ncpus=1,perfect=('mock' in dataset and perfect),extras=extras)
        _test.co_detect(colour,dc=dc,width=99,pp=False,ncpus=1,perfect=('mock' in dataset and perfect),extras=extras)
        print 'blah'
        

    if ncpus > 1:
        detection_list,job_server = pp_run(detections,run_function,pp_args,dependent_functions,modules,ncpus=ncpus,idle_delay=360.0,verbose=True,server=job_server)
    else:
        detection_list = []
        for d in detections:
            detection_list.append(d.co_detect(colour,dc=dc,width=99,pp=False,ncpus=1,perfect=('mock' in dataset and perfect),extras=extras))

    #this ^ runs co_detect with "form detections" set to True, so returns just the source detection files


    all_dets = []
    for s in detection_list:
        all_dets = all_dets+s

    args = [{'source[i]':None},{'perfect':False}]
    import detections as classes
    import detections
    run_function = 'get_input'
    dependent_functions = (selection_fn,slice_data,make_mask,make_detections,classes,getData,fits,detections)
    modules = ('clusters','numpy','parallel_qhull','operator','fits_table','operator','detections')
    detections,job_server = pp_run(list(all_dets),run_function,args,dependent_functions,modules,ncpus=ncpus,idle_delay=360.0,server=job_server,verbose=verbose)

    for d in detections:
        uber_det.src_data[d.colour_dna()] = d.cuts
        try:
            d.selection['overrides']['return_detection'] = True
        except:
            d.selection['overrides'] = {}
            d.selection['overrides']['return_detection'] = True
        if perfect:
            d.selection['overrides']['make_perfect'] = True

    run_function = 'parallel_detect'
    args = [{'source[i].cuts':None},{'source[i].selection':None}]
    dependent_functions = (deploy_vertices,qhull,parallel_detect,sel_edit,sci.richness,cherrypick,dataIO.Progress,sel_edit2)
    modules = ('science','parallel_detect','clusters','numpy','parallel_qhull','parallel_vertices','vperc','time','detections')
#    _ddd = parallel_detect(detections[0].cuts,detections[0].selection)
    detections,job_server = pp_run(detections,run_function,args,dependent_functions,modules,ncpus=ncpus,idle_delay=360.0,verbose=True,server=job_server)

    if permit_dropouts:
        cauldron = []
        all_clgals = {}
        all_second_stage_clusters = []
        #what is next colour?
        from operator import indexOf
        _colours = first_stage[0].selection['detection']['sequence']
        _ii = indexOf(_colours,colour)
        _ii += 1
        next_colour = _colours[_ii]
        next_magnitude = next_colour[-1]
        missing_in_detB = []
        for det in detections:     #detections is now the 2nd-stage detections :-/
            all_second_stage_clusters += list(det.clusters)
            if len(det.clusters) > 0:
                for g in det.cluster_galaxies:
                    all_clgals[g.extras['id']] = 1
        for galID in cl_gals.keys():
            try:
                dummy = all_clgals[galID]
            except KeyError:
                #this galaxy has been lost! It was in detA, but is not found in detB
                #does it have any nans for colours / magnitudes?
                #get me the galaxy object
                gal = galID_lookup[galID]
                if numpy.isnan(gal.magnitudes[next_magnitude]):
                    missing_in_detB.append(gal.extras['id'])
        missing_in_detB  = unique(missing_in_detB)
        check_these_clusters = []
        for galID in missing_in_detB:
            check_these_clusters = check_these_clusters + list(cl_gals_obj[galID])
        check_these_clusters = unique(check_these_clusters)  #because >1 missing_in_detB galaxy might belong to a cluster, so repeats
        if len(check_these_clusters) > 0:
            print "There are %d/%d clusters that need to be checked to ensure they are not 'impotent'"%(len(check_these_clusters),len(cl_members.keys()))
            descendent_lookup = {}
            for cluster in check_these_clusters:
                #get all galaxies in cluster
                descendent_lookup[cluster.dna] = []
                gals = cluster.galaxies
                invalid_mags = [numpy.isnan(gal.magnitudes[next_magnitude]) for gal in gals]
                first_stage_gal_ids = [g.extras['id'] for g in cluster.galaxies]
                nbad = countOf(invalid_mags,True)
                                #if Ngood <= 5 ---> cauldron

                if len(gals) - nbad < detections[0].selection['detection']['smallest_cluster']:
                    cauldron.append(cluster)
                else:
                    #look at all clusters in 2nd stage
                    #for each 2nd stage cluster, count no of galaxies from this 1st stage cluster

                    has_galaxies = []    #a 2nd stage cluster has X galaxies from 1st stage 
                    for cluster_stage2 in all_second_stage_clusters:
                        cs2_gal_ids = [g.extras['id'] for g in cluster_stage2.galaxies]
                        nmatch = 0
                        for id in first_stage_gal_ids:
                            if id in cs2_gal_ids:
                                nmatch += 1
                        if nmatch != 0:
                            has_galaxies.append(nmatch)
                            descendent_lookup[cluster.dna].append(cluster_stage2)
                    
                    if (len(has_galaxies) == 0) or (max(has_galaxies) < detections[0].selection['detection']['smallest_cluster']):
                        #1st_stage_gals in any 2nd_stage_cls < 5 ---> cauldron                        
                        #cluster (or "remnant of cluster")) was not identifed in 2nd stage detections
                        cauldron.append(cluster)  
                    else:
                        #this has been "detected" in the 2nd stage - don't add to cauldron
                        a = 1

                        
                        
        if len(cauldron) > 0:
            print 'There are %d/%d first-stage clusters that need to be included directly in 2nd-stage cherrypicking'%(len(cauldron),len(cl_members.keys()))
            a = 1

            #group all clusters in the cauldron by common colour intercept
            common_cm20 = {}
            cm20s = unique([cl.colour_dna() for cl in cauldron])
            for cm20 in cm20s:
                common_cm20[cm20] = []
            for cl in cauldron:
                common_cm20[cl.colour_dna()].append(cl)

            #now generate a set of detections, each with these cm20s, and add the clusters to them
            cauldron_detections = []
            for key in common_cm20.keys():
                det_tmp = common_cm20[key][0].parent.clone()
                det_tmp.addclusters(common_cm20[key])
                for cluster in det_tmp.clusters:
                    cluster.origin = 'cauldron'
                cauldron_detections.append(det_tmp)

            detections = detections + list(cauldron_detections)
            for d in cauldron_detections:
                if d.colour_dna() in uber_det.src_data.keys():
                    print 'Error - attempted to register cauldron detection, but already exists in uber_det.src_data (colour_dna=%d)'%(d.colour_dna())
                    raise SystemExit(0)
                uber_det.src_data[d.colour_dna()] = d.cuts

                

            #do stuff here - where count cl.dna > 1 in cauldron
            #dummy det?
            #det.add_clusters(cauldron)??
            #add it to some list
            

                
    #group the output 2nd stage detections by common intercept in the "base" (ie 1st) colour:
    c1 = unique([d.selection['colours'][colour]['intercept'] for d in detections])
    co_dets = {}
    for int in c1:
        co_dets[int] = []
    for d in detections:
        co_dets[d.selection['colours'][colour]['intercept']].append(d)
    ints = list(co_dets.keys())
    ints.sort()
    filters = []
    for i in ints:
        filters.append(co_dets[i])

    run_function = 'cherrypick'
    args = [{'source[i]':None},{'sci':sci},{'metric':detections[0].selection['detection']['metric']},{'cluster_zoo':False},{'include_solo':True},{'perfect':False}]
    all_filters = filters
    if perfect:
        pc_dummies = []
        for f in filters:
            sub_filters = []
            for sf in f:
                pc_det = det(None,sf.selection,dummy=True)
                pc_det.addclusters(sf.perfect_clusters,verbose=False)
                sf.perfect_clusters = []
                sf.perfect_cluster_galaxies = []
                try:
                    pc_det.selection['overrides']['as_perfect'] = True
                except:
                    pc_det.selection['overrides'] = {}
                    pc_det.selection['overrides']['as_perfect'] = True
                sub_filters.append(pc_det)
            pc_dummies.append(sub_filters)
        all_filters = list(filters)+list(pc_dummies)

#    collapsed_clusters = [cherrypick(slices,sci,metric=detections[0].selection['detection']['metric'],include_solo=True,perfect=False,verbose=False) for slices in all_filters]
    print 'Running cherrypicker..'
    

    collapsed_clusters_pp,job_server = pp_run(all_filters,run_function,args,dependent_functions,modules,ncpus=ncpus,idle_delay=36000000.0,verbose=True,server=job_server)

    junk = """\
    import os
    from dataIO import save
    if os.path.isfile('all_filters.dat') == False:
        save(all_filters,'all_filters.dat')
    """
    cp_stat = Progress(len(all_filters))
    collapsed_clusters = []
    for slices in all_filters:
        collapsed_clusters.append(cherrypick(slices,sci,metric=detections[0].selection['detection']['metric'],include_solo=True,perfect=False,verbose=False))
        cp_stat.update(0)

    true_clusters = []
    perfect_clusters = []
    t1 = iter(filters)
    t2 = iter(collapsed_clusters)
    while 1:
        try:
            f = t1.next()
            true_clusters.append(t2.next())
        except:
            break
    t1 = iter(filters)
    while 1:
        try:
            f = t1.next()
            perfect_clusters.append(t2.next())
        except:
            break

    if perfect == False:
        perfect_clusters = [False for p in filters]

    if perfect and (len(perfect_clusters) != len(true_clusters)):
        print 'mismatch'
        raise SystemExit(0)

    co_dets = []
    for i in xrange(len(filters)):
        _det = filters[i][0]
        colour_dna = unique([_cl.colour_dna() for _cl in true_clusters[i]])
        _det.src_data = {'total':uber_det.src_data['total']}
        for dna in colour_dna:
            _det.src_data[dna] = uber_det.src_data[dna]

        _det.addclusters(true_clusters[i],perfect=perfect_clusters[i],verbose=False)
        co_dets.append(_det)

    
    if richness:
        return co_dets,uber_det
#        co_dets = pp_richness_launcher(co_dets,uber_det,job_server=job_server)

    return co_dets
        
        
def pp_richness_launcher(dets,uber_det=None,job_server=None,construct=False,uber=False,ncpus=16,rfrac=0.8):
    from parallel_detect import parallel_detect
    from science import cherrypick,science
    import dataIO
    from dataIO import save,load,pp_run,selection_fn,slice_data,make_mask,make_detections,getData,sel_edit,sel_edit2,Progress
    import fits_table as fits
    from detections import detections as det

    if construct and not uber_det and len(dets) == 1:
        #this is a merged dataset - need to split
        cdna = {}
        cdna_sel = {}
        _det = dets[0]
        for cl in _det.clusters:
            cl_dna = cl.colour_dna()
            if cdna.has_key(cl_dna):
                cdna[cl_dna].append(cl)
            else:
                cdna[cl_dna] = [cl]
            cdna_sel[cl_dna] = cl.selection
        


    if construct or len(uber_det.src_data.keys()) <= 1:
        print 'Constructing datasources'
        if not uber_det.src_data.has_key('total'):
            uber_det.src_data = {}
            for c in uber_det.selection['detection']['sequence']:
                uber_det.selection['colours'][c]['width'] = 100.0
            uber_det.get_input()
            uber_det.src_data['total'] = uber_dets.cuts
        for d in dets:
            uber_det.src_data[d.colour_dna()] = d.cuts
            del d.cuts
        for d in dets:
            d.src_data = {'total':uber_det.src_data['total']}
            colour_dna = unique([_cl.colour_dna() for _cl in d.clusters])
            for dna in colour_dna:
                d.src_data[dna] = uber_det.src_data[dna]
            
    print 'Calculating cluster richness'
    bgcR,r80,bcg_centre,do_all = False,False,False,True
    args = [{'source[i]':None},{'bcgR':bgcR},{'r80':r80},{'bcg_centre':bcg_centre},{'do_all':do_all},{'uber':uber},{'rfrac':rfrac}]
    dependent_functions = ()
    modules = ('science','parallel_detect','clusters','numpy','parallel_qhull','parallel_vertices','vperc','time','detections')
    run_function = 'cluster_simple_richness'
#    oo = cluster_simple_richness(dets[0],do_all=True)
    completed,job_server = pp_run(list(dets),run_function,args,dependent_functions,modules,ncpus=ncpus,idle_delay=360.0,verbose=True,server=job_server)
    print 'Richness completed'
    return completed

