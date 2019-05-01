%% read CS FLASH kspace data
function [kSpace] = read_CSFlash_data(pathname)

%---------------
%pathname = ['/home/kostya/Documents/BCS_matlab_code/BCS_pros/48'];
method = readmethod(pathname);
acq = readacqp(pathname);

%---Number of Repetitions
NR = acq.NR;
%---Number of Slices
NSlices = acq.NSlices;
%---Number of Phase & RO steps
PE_steps = acq.size(2);
RO_steps = method.kSize(1);
%---Reduction Factor
RFactor = method.Rfactor;
%---Some SHIT-------
PE_total = PE_steps*RFactor;

%------Reading FID word_by_word
FID=fopen([pathname,'/fid']);
if(strcmp(acq.wordsize,'_32_BIT'));wSize = 'int32';else wSize = 'int16';end
RawData = fread(FID,inf,wSize,0,'l');
fclose(FID);

%-------Creating complex data from words-----
RawDataComplex = RawData(1:2:end) + sqrt(-1)*RawData(2:2:end);
[TotalPoints, ~] = size(RawDataComplex); 

%--Computing an ADC sampling size (num of points/words sampled at ones)
%--should be equal or bigger then RO and be
%--an integer value within ADC sampling range.
ADC = (0:1:31);
bit = 1;
while RO_steps > power(2, ADC(bit)); bit=bit+1; end

ADC_SampPoints = power(2, ADC(bit));

%--Computing k-space lines opposite permute procedure
%--FYI - sampling pocedure shuffles kLines but Gradp in ACQP contains 
%--relative position of the line by Index
index = method.Gradp2;
kLineIndex =round(PE_total*(index+1)/2 + 1);
kLineIndex=reshape(kLineIndex,[NSlices,PE_steps,NR]);
%kLineIndex = permute(kLineIndex,[2 3 1]);

if TotalPoints == ADC_SampPoints*NSlices*PE_steps*NR
    
    ADCMatrix = reshape(RawDataComplex,ADC_SampPoints,NSlices,PE_steps,NR);
    RawMatrix = ADCMatrix(1:RO_steps,:,:,:);

  %  [m2, kLineIndex] = mask_cs2d(pathname);
    
    kSpace = zeros(size(RawMatrix));
    for i=1:acq.NSlices
        for j=1:acq.NR
            kSpace(:,i,kLineIndex(i,:,j),j) = RawMatrix(:,i,:,j);
        end
    end
   % m2_1 = shiftdim(m2,-1);
  %  m2 = repmat(m2,[method.imSize(1),1,1]);
end