function status=odetpbar(t,y,flag, M0, Charge, Fields, Scaling)
    persistent tf tstart;
    
    if isempty(flag)
        % Integration steps
        ts=mean(abs(t));
        progress=100*ts/tf;
        textprogressbar(progress);
        status=0;
    else
        switch flag
            case 'init'     % Initializing progress bar
                tstart=tic;
                tf=max(abs(t));
                textprogressbar('ODE integration: ');
            case 'done'     % Finishing status function
                tf=[];
                textprogressbar('');
                display([ '   Integration time: ' num2str(toc(tstart))]);
                tstart=[];
            otherwise
                error('odetpbar:UnknownError',...
                    'Unknown error has occured');
        end
    end
end