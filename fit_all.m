fun = @(beta,x) (log(beta(1)+beta(2)*x(:,1).^beta(3)));
fun2 = @(beta,x) (beta(1)+beta(2)*x(:,1).^beta(3));
good_tracks = extract_good_tracks(alldata_control);
[lag_vec,msd,msd_err] = ensemble_msd_new(good_tracks,60,0.25);
figure(3);
clf;
errorbar(lag_vec,msd,msd_err,'b','linewidth',2);
set(gca,'xscale','log','yscale','log');
axis square;
hold on
[beta_control_all,R,J,covB,MSE]=nlinfit(lag_vec(1:60),log(msd(1:60,1)),fun,[0.003,0.001,0.5]);
plot(lag_vec(:,1),fun2(beta_control_all,lag_vec(:,1)),'k--','linewidth',2)
tracklen=cellfun(@length,good_tracks);
idx = find(tracklen>120);
tr{1} = good_tracks;
[fitpars, e_msd] = msd_resample(idx,tr);
fitpars
ddum={};
for ii=1:length(idx)
    ddum{ii} = good_tracks{idx(ii)};
end
[beta_boostrap_control] = construct_bootstrap(ddum);
[lag_vec,msd,msd_err] = ensemble_msd_new(ddum,60,0.25);
figure(4);clf;
errorbar(lag_vec,msd,msd_err,'b','linewidth',2);hold on
axis square;
set(gca,'xscale','log','yscale','log');
plot(lag_vec(1:30),fun_v(fitpars),'g-.','linewidth',2);
[beta_control,R,J,covB,MSE]=nlinfit(lag_vec(1:60),log(msd(1:60,1)),fun,[0.003,0.001,0.5]);
plot(lag_vec(:,1),fun2(beta_control,lag_vec(:,1)),'k--','linewidth',2)
%
good_tracks = extract_good_tracks(alldata_dznep);
[lag_vec,msd,msd_err] = ensemble_msd_new(good_tracks,60,0.25);
figure(3)
errorbar(lag_vec,msd,msd_err,'r','linewidth',2);
hold on;
[beta_dznep_all,R,J,covB,MSE]=nlinfit(lag_vec(1:60),log(msd(1:60,1)),fun,[0.003,0.001,0.5]);
plot(lag_vec(:,1),fun2(beta_dznep_all,lag_vec(:,1)),'m--','linewidth',2)
tracklen=cellfun(@length,good_tracks);
idx = find(tracklen>120);
tr{1} = good_tracks;
[fitpars, e_msd] = msd_resample(idx,tr);
fitpars
ddum={};
for ii=1:length(idx)
    ddum{ii} = good_tracks{idx(ii)};
end
[lag_vec,msd,msd_err] = ensemble_msd_new(ddum,60,0.25);
figure(4);
errorbar(lag_vec,msd,msd_err,'r','linewidth',2);
[beta_dznep,R,J,covB,MSE]=nlinfit(lag_vec(1:40),log(msd(1:40,1)),fun,[0.003,0.001,0.5]);
plot(lag_vec(:,1),fun2(beta_dznep,lag_vec(:,1)),'m--','linewidth',2);
plot(lag_vec(1:30),fun_v(fitpars),'c-.','linewidth',2);
[beta_boostrap_dznep] = construct_bootstrap(ddum);

function beta = construct_bootstrap(tracks)
    fun = @(beta,x) (log(beta(1)+beta(2)*x(:,1).^beta(3)));
    num_tracks = length(tracks);
    for iensemble=1:200
        idx = randperm(num_tracks);
        ddum = {};
        for ii = 1:round(0.8*num_tracks)          
            ddum{ii} = tracks{idx(ii)};
        end
        [lag_vec,msd,~] = ensemble_msd_new(ddum,60,0.25);
        beta(iensemble,:)=nlinfit(lag_vec(1:40),log(msd(1:40,1)),fun,[0.003,0.001,0.5]);
    end
end

