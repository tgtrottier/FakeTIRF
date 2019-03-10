%% File parameters
fileDirectory = cat(2,pwd,'\Generated Data\');
fileName	  = 'BleachVirialExpVaryBigMass3r1-3'; %Omit extension

%% Render parameters
renderPreview    = true;  %Whether to show user the plots as they are generated
renderXPixels    = 512;   %X pixel count
renderYPixels    = 512;   %Y pixel count
renderResolution = .27e-6;  %Image resolutuon (meters)
renderTimeStep   = 1/32; %Time step between frames (seconds)
renderDuration   = 10;     %How long to render for (seconds)
renderOversize   = 2;     %How much larger render for particle positions to account from them walking in and out of view

%% Particle physical properties
density = 2633; % density in kg/m^3
particleCount      = 256;   %How many to draw (may be increased based on renderOversize)
particleBrightness = 1;  %The maximum brightness of each particle [0,1]. 1 saturates image. Particles on top of each other add intensity.
bleachConst = 0.05; %time constant for particle bleaching
%particleRadius     = .5e-6; %The radius of the spot in the image (meters) (spots standard deviation)
particleRadius = (1+1.*(randn(particleCount*renderOversize^2,1)))*1e-7; % create a vector of radii with (avg+stddev)*scale 
bigBoyFactor  = 0.00;%generate a random array of BIG boys and add them to the particle list [0,1]
particleRadius = max((5+1.*(randn(particleCount*renderOversize^2,1)))*1e-6.*(rand(particleCount*renderOversize^2,1)<bigBoyFactor),particleRadius);
% particle speed is equal to the average + a normally distributed portion
% around the average with a variance of 9kbt/4pi*rho*r^3 consistent with
% the virial theorem.
virialScale = 1/800; % virial scaling factor required for unaccounted forces
particleSpeed  =  virialScale*particleRadius.^(-3/2)*sqrt(9*1.38e-23*298/(4*pi*density))+particleRadius.^(-3/2)*9*1.38e-23*298/(4*pi*density).*randn(particleCount*renderOversize^2,1);%Average speed of particles (meters per second) (mean squared speed)

%% Extraneous
%Define how 1/e^2 distance for brightness for the particles. This should
%be the same distance at which the electric field decays to 1/e.
%Set to inf to disable the brightness dimming.
physicsBrightnessFalloffRate = 2400e-9 ; %(meters)
%This is the decay function for the brightness, which uses the previously
%defined parameter. This function probably won't need to be altered.
physicsBrightnessFalloffFunc = @(distances) exp(-2*distances/physicsBrightnessFalloffRate);
%Detachment rate constant
%parameter to decide rate constant of detachment
%The parameter is multiplied by the time step and detachesa based on the
%function rand < detachConst*renderTimeStep
detachConst = 3;
%sticking parameter number of radii inside which the particle will get stuck to
%the surface
stickP = 1;
%noise intensity and variance [0,1]
noiseI=0.01;noiseV=0.01;





%% Perform file save prep
%Make directory if needed, but prevent warnings if the folder is already
%present.
[success,~,~] = mkdir(fileDirectory);
if ~success
	error('Could not create/access the desired animations save directory.')
end; clear success
fullFilePathTiff = cat(2,fileDirectory,fileName,'.tiff');
fullFilePathMat  = cat(2,fileDirectory,fileName,'.mat');

%If the file already exists, ask the user if they want to overwrite it
if numel(dir(fullFilePathTiff)) == 1 || numel(dir(fullFilePathMat)) == 1 %file already exists

	response = questdlg('Previously saved files already exists. Overwrite them?', ...
		'Filename Collision',...
		'Yes, overwrite.','No, halt function.','No, halt function.');
	if ~strcmp(response,'Yes, overwrite.')
		error('Function stopping because of name collision. Change fileName parameter.')
	end
end; clear response

%% Perform render prep
%Positions which are sampled in image
posX =      ((0:renderXPixels-1) - (renderXPixels-1)/2) * renderResolution;
posY = flip(((0:renderYPixels-1) - (renderYPixels-1)/2) * renderResolution); %Flip Y so it plots correctly

%Plan out which times will be measured
times = 0 : renderTimeStep : renderDuration;
numFrames = numel(times);

%Define a function which will be used to bound, and another which will
%handle the scaling and rounding of intensity values to match the desired
%index range.
unitBound = @(x) max(min(x,1),0);
brightnessToIndex = @(x) round(unitBound(x)*255);

%% Evaluate where all the particles are
[centerX,centerY,centerZ] = GeneratePositions(particleCount,renderOversize,posX,posY,numFrames,renderTimeStep,particleSpeed,detachConst,stickP,particleRadius);
save(fullFilePathMat,'centerX','centerY','centerZ','posX','posY')

%% Render the frames
%Create figure, but only if the user wants a preview
if renderPreview
	close all;
	figure();
	%Specify the image should be draw with the gray colormap
	colormap('gray')
end

%Loop over the different frames to render
for timeInd = 1:numFrames
	%Calculate the intensities for the frame. Bound the intensities to
	%[0,1] before the addition of noise so the noise is consistently
	%applied.
	intensities = imnoise(...
		unitBound(...
			RenderFrame(...
				posX,...
				posY,...
				centerX(:,timeInd),...
				centerY(:,timeInd),...
				particleRadius,...
				particleBrightness*exp(-bleachConst*timeInd*renderTimeStep)*physicsBrightnessFalloffFunc(centerZ(:,timeInd))...
			)...
		),...
		'gaussian',noiseI,noiseV);
	%If the user wants to see previews,
	if renderPreview
		%Draw that frame with gaussian noise. The brightness saturation is
		%set to 1.
		imagesc(intensities,[0,1]);
		daspect([1,1,1]); %Lock aspect ratio to 1:1
		drawnow;
		pause(1/30);
	end
	%If its the first frame, save to the file anew
	if timeInd == 1
		imwrite(brightnessToIndex(intensities),gray(256),fullFilePathTiff,'WriteMode','overwrite');
	%Otherwise append the frame
	else
		%I have noticed that this process can be so fast that the file
		%write doesnt have time to finish. Im adding this so that it will
		%catch if it fails and will try again until succeeds.
		breakCount = 1;
		while true
			try
				imwrite(brightnessToIndex(intensities),gray(256),fullFilePathTiff,'WriteMode','append');
				break;
			catch
				%All of these counters and fprints are just to help the
				%user know how many times it failed in a row. This also
				%implements a way to delete some previously written
				%characters in the console so the message is concise even
				%if it runs hundreds of times.
				if breakCount == 1
					fprintf('The tiff frame write was too fast (frame %u,  ',timeInd);
					charCount = 1; %Added a space to the last string, better to remove one character than zero...
				end
				breakCount = breakCount + 1; %Do this first since 1 failure means 2 attempts (1 successful), etc.
				fprintf(repmat('\b',1,charCount)); %Remove characters as needed
				charCount = fprintf('%u attempts)\n',breakCount);
			end
		end
	end
end




