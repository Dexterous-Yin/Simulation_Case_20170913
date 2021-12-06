%{
This file is used to plot the energy spectrum observed by RBSP, and
simulate the observation using the result from Trace_magdip_proton_rbsp_case.m
Zefan -- 21-03-18
%}

%% plot the RBSP observation
filename = 'data\rbspb_rel04_ect-hope-pa-l3_20170912_v7.4.0';
info = cdfinfo(filename);
out = cdfread(filename);
out = out(1:3434,:);
filename_add = 'data\rbspb_rel04_ect-hope-pa-l3_20170913_v7.2.0';
info_add = cdfinfo(filename_add);
out_add = cdfread(filename_add);
out_add = out_add(1:3498,:);
%
% 10-Epoch-ion % 47-FPDU % 28-PA %16-Energy_Ion %6-Energy_delta
time = zeros(1,length(out)+length(out_add));
hope_energy = out{1,16}./1000; %[keV]
hope_energy_delta = out{1,6}./1000; %[keV]
xE = 1:72;
xE_interp = 0.5:1:72.5;
p = polyfit(xE,log(hope_energy),1);
hope_energy_interp = exp(polyval(p,xE_interp));
%
Estart = 68;
Eend = 71;
Enum = Eend-Estart+1;
PAid = find(out{1,28}==90);
PAid_sPA = find(out{1,28}==18 | out{1,28}==162 | out{1,28}==36 | out{1,28}==144);
hope_flux = zeros(length(out)+length(out_add),Enum+1);
hope_flux_sPA = zeros(length(out)+length(out_add),Enum+1);
hope_E = hope_energy_interp(Eend+1:-1:Estart)'; % from High-E to Low-E
%
for i=1:length(out)
    for j=1:Enum+1
        time(i)=todatenum(out{i,10});
        hope_flux(i,j) = out{i,47}(PAid,Eend-j+1);
        temp_flux_spA = out{i,47}(PAid_sPA,Eend-j+1);
        temp_flux_spA(temp_flux_spA<=0)=nan;
        hope_flux_sPA(i,j) = mean(temp_flux_spA,'omitnan');
    end
end
for i=1:length(out_add)
    for j=1:Enum+1
        time(i+length(out))=todatenum(out_add{i,10});
        hope_flux(i+length(out),j) = out_add{i,47}(PAid,Eend-j+1);
        temp_flux_spA = out_add{i,47}(PAid_sPA,Eend-j+1);
        temp_flux_spA(temp_flux_spA<=0)=nan;
        hope_flux_sPA(i+length(out),j) = mean(temp_flux_spA,'omitnan');
    end
end

%%
% clear out;
formatIn = 'yyyy-mm-dd/HH:MM:SS';
start_time = datenum('2017-09-12/23:30:00',formatIn);
end_time = datenum('2017-09-13/03:30:00',formatIn);
trange = find(time>=start_time & time<=end_time);
trange=[trange,trange(end)+1];
time_pco = time(trange);
flux_sim = log10(hope_flux(trange,:)');
flux_sim_sPA = log10(hope_flux_sPA(trange,:)');

figure('Color',[1 1 1]);
subplot(4,1,1);
pco_hope = pcolor(time_pco,hope_E,flux_sim);
set(pco_hope,'LineStyle','none');
set(gca,'YScale','log');
% set(gca,'YDir','reverse');
colormap jet;
c_hope = colorbar;
c_hope.Label.String = {'Log10 Flux';'90PA'};
caxis([4,6]); 
xlim([start_time,end_time]);
ttick = start_time:datenum('0000-01-00/00:30:00',formatIn):end_time;
tticklabel = datestr(ttick,'HHMM');
set(gca,'XTick',ttick);
set(gca,'XTickLabel',tticklabel);
ylabel({'Energy';'keV'});
title('Observation');
set(gca,'Layer','top');

subplot(4,1,2);
pco_hope_sPA = pcolor(time_pco,hope_E,flux_sim_sPA);
set(pco_hope_sPA,'LineStyle','none');
set(gca,'YScale','log');
colormap jet;
c_hope_sPA = colorbar;
c_hope_sPA.Label.String = {'Log10 Flux';'small PA'};
caxis([4,6]); 
xlim([start_time,end_time]);
% xlim([datenum('2019-02-27/12:00:00',formatIn),datenum('2019-02-27/16:00:00',formatIn);]);
set(gca,'XTick',ttick);
set(gca,'XTickLabel',tticklabel);
ylabel({'Energy';'keV'});
set(gca,'Layer','top');
%%
% plot simulation results
time_source = (22+47/60)*60*60;
load('simulation_result\sim_90_X5.5_Lm_full.mat');
simflux = zeros(length(E_ob),length(time_ob));
for i=3:6
    for j=1:length(time_ob)
        whichm = 3; % AP8max
        energy_bk = E_ob(i)*1e-3;
        BBo = 1;
        L_bk = L_ob(j);
        bkground = onera_desp_lib_get_ae8_ap8_flux(whichm,energy_bk,BBo,L_bk)/1000/4/pi;
        inject = 0;
        if flux_ob(i,j)>=time_source
            tempE = W_rec{i,j}(end); % energy at source region
            source = 3e5*(10/tempE)^0;  % flux at source
            inject = source/tempE*W_rec{i,j}(1);
        end
        simflux(i,j) = bkground+inject;
    end
end

% LOW resolution
% plot (NOTE the pco direciton has been reversed)
subplot(4,1,3);

xE = 1:length(E_ob);
xE_interp = 0.5:1:length(E_ob)+0.5;
p = polyfit(xE,log(E_ob),1);
E_ob_interp = exp(polyval(p,xE_interp));

E_sim = E_ob_interp(7:-1:3);
nE = length(E_ob);
flux_sim = log10(simflux(6:-1:2,:));
pco_sim = pcolor(time_ob,E_sim,flux_sim);
set(pco_sim,'LineStyle','none');
set(gca,'YScale','log');
% set(gca,'YDir','reverse');
colormap jet;
c = colorbar;
c.Label.String = {'Log10 Flux';'90PA'};
caxis([4,6]); 
% ylim([10,100])
xstart = 23.5;
xend = 27.5;
simttick = xstart:0.5:xend;
xlim([xstart,xend]);
formatIn = 'yyyy-mm-dd/HH:MM:SS';
tticklabel = datestr(datenum('2017-09-12/00:00:00',formatIn)+simttick./24,'HHMM');
set(gca,'XTickLabel',tticklabel);
title('Simulation');
ylabel({'Energy';'keV'});
set(gca,'Layer','top');


%%
% plot simulation results with low pitch angle
load('simulation_result\sim_27_X5.5_Lm_full.mat');
simflux = zeros(length(E_ob),length(time_ob));
for i=3:6
    for j=1:length(time_ob)
        whichm = 3; % AP8max
        energy_bk = E_ob(i)*1e-3;
        BBo = 1;
        L_bk = L_ob(j);
        bkground = onera_desp_lib_get_ae8_ap8_flux(whichm,energy_bk,BBo,L_bk)/1000/4/pi;
        inject = 0;
        if flux_ob(i,j)>= time_source
            tempE = W_rec{i,j}(end); % energy at source region
            source = 3e5; %exp(interp1(log(th_E),log(th_flux),log(tempE),'linear','extrap')); % flux at source
            inject = source/tempE*W_rec{i,j}(1);
        end
        simflux(i,j) = bkground+inject;
    end
end

% plot (NOTE the pco direciton has been reversed)
subplot(4,1,4);
xE = 1:length(E_ob);
xE_interp = 0.5:1:length(E_ob)+0.5;
p = polyfit(xE,log(E_ob),1);
E_ob_interp = exp(polyval(p,xE_interp));

E_sim = E_ob_interp(7:-1:3);
nE = length(E_ob);
flux_sim = log10(simflux(6:-1:2,:));
pco_sim = pcolor(time_ob,E_sim,flux_sim);
set(pco_sim,'LineStyle','none');
set(gca,'YScale','log');
% set(gca,'YDir','reverse');
colormap jet;
c = colorbar;
c.Label.String = {'Log10 Flux';'small PA'};
caxis([4,6]); 
% ylim([10,100])
xstart = 23.5;
xend = 27.5;
simttick = xstart:0.5:xend;
xlim([xstart,xend]);
formatIn = 'yyyy-mm-dd/HH:MM:SS';
tticklabel = datestr(datenum('2017-09-12/00:00:00',formatIn)+simttick./24,'HHMM');
set(gca,'XTick',simttick);
set(gca,'XTickLabel',tticklabel);
ylabel({'Energy';'keV'});
set(gca,'Layer','top');



