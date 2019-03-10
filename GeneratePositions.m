%This function generates the x,y,z positions of all the particles being
%modeled.
function [centerX,centerY,centerZ] = GeneratePositions(particleCount,renderOversize,posX,posY,numFrames,renderTimeStep,particleSpeed,detachConst,stickP,particleRadius)
	
	%Provide the same particle density that the user asked for, on average,
	%but increase the number of particles to be traced so that extra may be
	%simulated outside of the rendering window. This allows them to walk in
	%and out of frame without it being apparent that we are simulating a
	%finite domain.
	numParticles = ceil(particleCount * renderOversize^2);
	
	%Scale up the domain (used only for the initial positions) to the
	%proper size. posX, posY are assumed sorted.
	xDomain = mean(posX([1,end])) + [-1,1]/2*renderOversize*abs(diff(posX([1,end])));
	yDomain = mean(posY([1,end])) + [-1,1]/2*renderOversize*abs(diff(posY([1,end])));
	
	%Initialize the position storage containers. Each row represents a
	%different particle, and each column represents a different time (frame).
	centerX = nan(numParticles,numFrames);
	centerY = nan(numParticles,numFrames);
	centerZ = nan(numParticles,numFrames);
	
%FROM HERE ON IS REASONABLE TO EDIT
	
	%Populate the particles with initial conditions.
	centerX(:,1) = rand(numParticles,1).*diff(xDomain) + min(xDomain); %Uniformly place over xDomain
	centerY(:,1) = rand(numParticles,1).*diff(yDomain) + min(yDomain); %Uniformly place over yDomain
	centerZ(:,1) = 0; %Temporarily, place all at boundary. If this is changed, then we probably need more particles from the start
	
	%Loop over each particle and frame. Looping over both allows for
	%boundary constraints to be enforced at each frame.
	for frameInd = 2:numFrames %First frame is already defined by initial positions
		for particleInd = 1:numParticles
           
            %test to see if we can get particles stuck
            if centerZ(particleInd,frameInd-1) < renderTimeStep*stickP*particleRadius
                centerZ(particleInd,frameInd) = 0;
                centerX(particleInd,frameInd) = centerX(particleInd,frameInd-1) + particleSpeed(particleInd)/sqrt(2)*renderTimeStep*randn(1); %speed^2/3 is the variance of displacements in a single direction
                centerY(particleInd,frameInd) = centerY(particleInd,frameInd-1) + particleSpeed(particleInd)/sqrt(2)*renderTimeStep*randn(1);
                stuckThisF = true;
            end
                %check to see if you pop 
            if (centerZ(particleInd,frameInd-1) == 0 && rand < exp(-detachConst*renderTimeStep)) 
                centerZ(particleInd,frameInd) = centerZ(particleInd,frameInd-1) + particleSpeed(particleInd)/sqrt(3)*renderTimeStep*randn(1);
                centerX(particleInd,frameInd) = centerX(particleInd,frameInd-1) + particleSpeed(particleInd)/sqrt(3)*renderTimeStep*randn(1); %speed^2/3 is the variance of displacements in a single direction
                centerY(particleInd,frameInd) = centerY(particleInd,frameInd-1) + particleSpeed(particleInd)/sqrt(3)*renderTimeStep*randn(1);
                %if you dont pop off then keep moving in 2D
            elseif centerZ(particleInd,frameInd-1) == 0 
                centerZ(particleInd,frameInd) = 0;
                centerX(particleInd,frameInd) = centerX(particleInd,frameInd-1) + particleSpeed(particleInd)/sqrt(2)*renderTimeStep*randn(1); %speed^2/3 is the variance of displacements in a single direction
                centerY(particleInd,frameInd) = centerY(particleInd,frameInd-1) + particleSpeed(particleInd)/sqrt(2)*renderTimeStep*randn(1);
            else
                %if you arent stuck then move in 3d
                centerZ(particleInd,frameInd) = centerZ(particleInd,frameInd-1) + particleSpeed(particleInd)/sqrt(3)*renderTimeStep*randn(1);
                centerX(particleInd,frameInd) = centerX(particleInd,frameInd-1) + particleSpeed(particleInd)/sqrt(3)*renderTimeStep*randn(1); %speed^2/3 is the variance of displacements in a single direction
                centerY(particleInd,frameInd) = centerY(particleInd,frameInd-1) + particleSpeed(particleInd)/sqrt(3)*renderTimeStep*randn(1);
            end
            
            
          
            
            
            
            
		end
		%Enforce boundary constraints
		centerZ(:,frameInd) = max(centerZ(:,frameInd),0); %Enforce z>0
	end
end