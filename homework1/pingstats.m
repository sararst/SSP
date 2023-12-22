function out=pingstats(machine,num,verbose)

% function out=pingstats(machine,num,verbose)
%
% Example : stats=pingstats('isl.stanford.edu',100,'v')
%   or      stats=pingstats('www.google.com',100,'v')
%   
%   machine : name or IP address of the computer to be "pinged"
%   num     : number of ping operations
%   verbose : if different from '' (empty string) the ping state gets displayed
%
%  For Homework 1 of SSP.
%  Executes pings to computer "machine" "num" times.
%  Lost packets are ignored (as if no ping was launched).
%  The vector "out" contains the ping times in milliseconds.
%  L.D.


i=1;
if (ispc)
    str1='ping ';
    str2=machine;
    str3=' -n 1';    
else
    str1='ping ';
    str2=machine;
    str3=' -c 1 -s 64 -i 1 ';
end

str=[str1 str2 str3];
format compact;

while i<=num      
    [s,w]=unix(str);
    if ~isempty(strfind(w,'time='))
        %extract the value given in "time=xx ms"
        temp1 = regexp(w,'time=','split');
        temp2 = regexp(temp1{2},'ms','split');
        w     = char(temp2{1});
        
        ww=str2num(w);
        if ~(isempty(ww))
            out(i)=ww;
            if ~(isempty(verbose))
                disp(i);
            end;
            i=i+1;
        end;     
    end
end;

format; 
