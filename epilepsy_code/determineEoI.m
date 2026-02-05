function[EoI] = determineEoI(channel,time, channel_labels)
    %Compute the envelope
    y = hilbert(channel);
    env = abs(y);

    %Standard Deviation and Mean of the envelope
    env_mean = mean(env);
    env_SD = std(env);

    %Set the threshold
    thr = env_mean + 1* env_SD;
    duration_thr = 0.5 * thr;


    %Mark The EoI
    EoI = zeros(size(env));
    for i=1:length(env)
        if env(i)>=thr
            EoI(i) = channel(i);
        else
            EoI(i) = 0;
        end
    end

    %I dont want to study the EoI above the spike
    for i = 1:length(EoI)
        if time(i)>0.1
            EoI(i) = 0;
        end
    end


    %Look for the duration of the EoI
    for i=1:length(EoI)
        if EoI(i)~=0
            %fprintf('\nEoI in sample %d\n',i);
            value = 0;
            rest = 1;
            while value==0 && (i - rest)>=1 
                if env(i - rest)>= duration_thr || env(i - rest) <= -duration_thr 
                    EoI(i - rest) = channel(i - rest);
                    %fprintf(' we have add the sample %d\n',(i-rest));
                    rest = rest + 1;
                else
                    value = 1;
                end

            end
        end
    end

    for i=1:length(EoI)
        if EoI(i)~=0
            value = 0;
            sum = 1;
            while value==0 && (i + sum)<=length(EoI)
                if env(i + sum)>= duration_thr|| env(i + sum) <= -duration_thr 
                    EoI(i + sum) = channel(i + sum);
                    sum = sum + 1;
                else
                    value = 1;
                end
            end
        end
    end

    %Merge together EoI in an interval less than 30ms (10 samples in
    %between) (to be accurate 7.68 samples are 30 ms with a fs = 256)
        interval_threshold = (abs(time(1)) - abs(time(10)));
        i = 1;
        while i <= length(EoI)
            if EoI(i)~=0 && i<length(EoI) 
                a = 0;
                %fprintf('\ni before while loop:%d\n',i);
                z = i+1;
                cont = i;
                while a == 0 && z<=length(EoI)
                    if EoI(z)~=0
                        z = z + 1;
                        a = 0;
                    else
                        a=1;
                        i = z-1;
                        cont = z-1;
                    end
                end
                %fprintf('\ni after while loop:%d\n',i);
                %fprintf('cont:%d\n',cont);
                b = cont+1;
                c = 0;
                    while b<=length(EoI) && c==0
                        if EoI(b)~=0 
                            if time(b)<0 && time(cont)<0
                                interval = abs(time(cont)) - abs(time(b));
                            elseif time(cont)<0 && time(b)==0
                                interval = abs(time(cont));
                            elseif time(cont)<0 && time(b)>0
                                interval = abs(time(cont)) + time(b);
                            elseif time(cont)>0 && time(b)>0
                                interval = time(b) - time(cont);
                            else
                                interval = 10000;
                            end
                            %fprintf('For pos %d and pos %d, the interval is: %d\n',cont,b,interval);
                            c = 1;
                            if interval<=interval_threshold
                                for j = (cont+1):b
                                    if EoI(j) == 0
                                        EoI (j) = channel(j);
                                        %fprintf('we have include the EoI in pos %d because the interval is %d\n',j,interval);
                                    end
                                end
                            end
                        end
                        b = b + 1;
                    end
            end
            i = i+1;
        end
end
