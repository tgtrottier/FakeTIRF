%This function intelligently adds the gaussian spots to the image
function intensities = RenderFrame(posX,posY,centerX,centerY,radius,brightness)
	windowOversize = 5; %How many radius distance to draw out
	
	%Determine the number of spots that will need to be drawn.
	numSpots = numel(centerX);
	
	%If the user chooses to pass a scalar for radius or brightness, turn
	%them into the properly sized vector by repeating the value they passed
	if numel(radius) == 1
		radius = repmat(radius,1,numSpots);
	end
	if numel(brightness) == 1
		brightness = repmat(brightness,1,numSpots);
	end
	
	%Create a place to store the new image.
	intensities = zeros(numel(posY),numel(posX)); %Flipped x y so it matches imagesc format
	
	%Loop over the different bright spots.
	for spotInd = 1:numSpots
		%Find which pixels are within a rectangular neighborhood of the
		%center point
		keepX = abs(posX - centerX(spotInd)) <= radius(spotInd) * windowOversize;
		keepY = abs(posY - centerY(spotInd)) <= radius(spotInd) * windowOversize;
		
		%Create a small meshgrid of the positions corresponding to these
		%pixels.
		[subPosX,subPosY] = meshgrid(posX(keepX),posY(keepY));%meshgrid automatically flips the order to match the imagesc format
		
		%Add onto the existing intensity matrix by a gaussian with
		%corresponding center. This makes use of the subPixX and subPixY to
		%only write to a small portion of the points in the image. This
		%will drastically cut down on time spent to make each frame.
		if sum(keepX)*sum(keepY) ~= 0
			intensities(keepY,keepX) = ...
				intensities(keepY,keepX) ...
				 + ...
				brightness(spotInd) ...
				 * ...
				exp(...
					-(...
						(subPosX - centerX(spotInd)).^2 ...
						 + ...
						(subPosY - centerY(spotInd)).^2 ...
					)/(2*radius(spotInd)^2)...
				);
		end
	end
end

