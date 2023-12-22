%w=window(type,L);
%
%generates a symmetric window w (column vector) of length L
%the centerpoint of the window has height 1.
%type is 'bartlett' or 'blackman' or 'boxcar'(=rectangular) or 'hamming'
%	  or 'hanning' or 'triang'

function w=window(type,L)

w=ones(L,1);

if strcmp(type,'bartlett') 
	w=2*(0:(L-1)/2)/(L-1);
	if rem(L,2)
        	% It's an odd length sequence
        	w = [w w((L-1)/2:-1:1)]';
	else
        	% It's even
        	w = [w w(L/2:-1:1)]';
	end
  
elseif strcmp(type,'blackman')==1
	w=(.42 - .5*cos(2*pi*(0:L-1)/(L-1)) + .08*cos(4*pi*(0:L-1)/(L-1)))';
elseif strcmp(type,'boxcar')==1 w=ones(L,1);
elseif strcmp(type,'hamming')==1 w=.54 - .46*cos(2*pi*(0:L-1)'/(L-1));
elseif strcmp(type,'hanning')==1 w=.5*(1 - cos(2*pi*(1:L)'/(L+1)));
elseif strcmp(type,'triang')==1
	if rem(L,2)
        	% It's an odd length sequence
        	w = 2*(1:(L+1)/2)/(L+1);
        	w = [w w((L-1)/2:-1:1)]';
	else
        	% It's even
        	w = (2*(1:(L+1)/2)-1)/L;
        	w = [w w(L/2:-1:1)]';
	end

end;

