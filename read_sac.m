function [seismograph] = read_sac(sacfile,endin)

%%This function reads a sac file and turns it into a 3 column matrix
%%INPUT: sac file name
%%OUTPUT: matrix with the first column time, the second column the
%%amplitude and third column sac parameters



%%Here we open up the sac file
seismograph=[];
try
    if(strcmp(endin,'b'))
    fid = fopen(sacfile,'r','ieee-be');
    else
       fid = fopen(sacfile,'r','ieee-le');
    end
catch
    seismograph=[];
end

if(fid~=-1);
%%Read in header information
header(1:70) = fread(fid,70,'single');
header(71:110) = fread(fid,40,'int32');
header(111:302) = (fread(fid,192,'char'))';

%%Allocate seismograph matrix
seismograph=zeros(header(80),3);

%%Set the amplitudes
seismograph(:,2)=fread(fid,'single');

%%Input the time times
seismograph(:,1) = linspace(header(6),header(6)+header(7),header(80)); 


%%Here we store all of the header information
%%Here is the sample rate
seismograph(1,3)=header(1);
%%Arrival marker
seismograph(2,3)=header(9);

%%Station information
%%Station latitude
seismograph(3,3)=header(32);
%%Station longitude
seismograph(4,3)=header(33);
%%Station elevation
seismograph(5,3)=header(34);
%%Station depth
seismograph(6,3)=header(35);

%%Event information
%%Event latitude
seismograph(7,3)=header(36);
%%Event longitude
seismograph(8,3)=header(37);
%%Event depth
seismograph(9,3)=header(39);
%%Event magnitude
seismograph(10,3)=header(40);
%%Event distance
seismograph(11,3)=header(51);
%%Event azimuth
seismograph(12,3)=header(52);
%%Event great circle arc
seismograph(13,3)=header(54);
%%Event year
seismograph(14,3)=header(71);
%%Event  day
seismograph(15,3)=header(72);
%%Event hour
seismograph(16,3)=header(73);
%%Event minute
seismograph(17,3)=header(74);
%%Event second
seismograph(18,3)=header(75);
%%Event second fraction
seismograph(19,3)=header(76);

seismograph(20,3)=header(12);


fclose(fid);
end

