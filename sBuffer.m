function [varargout] = sBuffer(func, varargin)
%SBUFFER Buffers output value of a function for speed up.
%If all input arguments match a previous result, return that result
%instead of calling the function itself.

persistent tBufferIn;
persistent tBufferOut;

if(nargout==0)
   nargouts=1; 
else
  nargouts=nargout;
end

if isempty(tBufferIn) %first call
    varargout{1:nargouts}=feval(func,varargin{:});
    tBufferIn= varargin;
    tBufferOut= varargout;
end

for i=1:nargin-1
    if isempty(varargin{i})
        varargout={NaN};
        return;
    elseif  ischar(varargin{i})
      indx(:,i)=strcmp(tBufferIn(:,i),varargin{i});  
    else
      indx(:,i)=([tBufferIn{:,i}] == varargin{i});   
    end  
end


%if one row meets all conditions
   if( any(all(indx')) )
       varargout = tBufferOut(all(indx'),:);
   else
       varargout{1:nargouts}=feval(func,varargin{:});
       
       
       if length(varargin{1})==1 %only save if identificator is unique
           tBufferIn = [tBufferIn; varargin];
           tBufferOut = [tBufferOut; varargout];
       end
       
   end


end

