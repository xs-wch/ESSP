function mycnt = myEEGfilter(filename)
close all
cntfile = [filename,'.cnt'];
mycnt=loadcnt(cntfile,'dataformat','int32');
srate=mycnt.header.rate;
mycnt.data = eegfilt(mycnt.data,srate,1,0,0,5*srate);
mycnt.data = eegfilt(mycnt.data,srate,0,40,0,srate);
savefile = [filename,'_filterd'];
save(savefile, 'mycnt');
