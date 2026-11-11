using LinearAlgebra
using StaticArrays
using Random
using Plots


macro swap(x,y)
    quote
       local tmp = $(esc(x))
       $(esc(x)) = $(esc(y))
       $(esc(y)) = tmp
    end
end

@inbounds function avalanche!(origin :: Int64, activeSites, activeSitesNew, isActiveNew, fth, f, g, k0, k1, k2 )
    L = size(fth)[1]
    nActiveSites = 1;
    activeSites[1] = origin;
    S_tot = 0;
    while nActiveSites > 0
        S_tot  = S_tot + nActiveSites;
        for k in 1:nActiveSites 
            site = activeSites[k];
            g[site] -= (2*k1+k0);
            f[site] -=  2*k2;
            fth[site] = randexp()^(1/gamma_rnd);
            for shift in (-1,1)
                nn = site + shift;
                if nn == 0
                    nn = L;
                end
                if nn == L+1
                    nn = 1;
                end
                g[nn] += k1;
                f[nn] += k2;
            end
        end
        nActiveSitesNew = 0;
        for k in 1:nActiveSites 
            site = activeSites[k];
            for shift in (-1,0,1)
                check_site = site + shift
                if check_site == 0
                    check_site = L
                elseif check_site == L+1
                    check_site = 1
                end
                delta = fth[check_site] - f[check_site] - g[check_site];
                if (delta < 0)
                    if (isActiveNew[check_site] == false)
                        isActiveNew[check_site] = true;
                        activeSitesNew[nActiveSitesNew+1] = check_site;
                        nActiveSitesNew += 1;
                    end
                end
            end 
            activeSites[k] = 0;
        end
        for k in 1:nActiveSitesNew
            site_to_reset = activeSitesNew[k];
            isActiveNew[site_to_reset] = false;
        end
        @swap(activeSites, activeSitesNew);
        nActiveSites = nActiveSitesNew;
    end
    return S_tot
end


function trigger!(fth, f, g)
    L = size(thresholds)[1];
    min_drive = Inf;
    origin_drive = 0;
    max_relax = 0.0;
    origin_relax = 0;
    for i in 1:L
        num = fth[i]-g[i];
        den = f[i];
        if (num < 0) && (den < 0)
            relax_amount = num/den;
            if (relax_amount < 1.0) && (relax_amount > 0) && (relax_amount > max_relax)
                max_relax = relax_amount;
                origin_relax = i;
            end
        end 
        drive_amount = fth[i] - g[i];
        if(drive_amount < min_drive)
            min_drive = drive_amount;
            origin_drive = i;
        end  
    end
    if (origin_relax != 0)
        for i in 1:L
            f[i] = f[i]*max_relax;
        end
        f[origin_relax] = fth[origin_relax] - g[origin_relax];
        return origin_relax, true, max_relax;
    else
        for i in 1:L
            f[i] = 0.0;
            g[i] += min_drive;
        end 
        g[origin_drive]  = fth[origin_drive];
        return origin_drive, false, min_drive; 
    end 
end

seed = 1101252;
Random.seed!(seed);

L=1024;
gamma_rnd = 2.0
k1 = 0.1;
k2 = 0.1;
k0 = 0.001;
thresholds = randexp(Float64,L);
thresholds = thresholds.^(1.0/gamma_rnd);
viscoel = zeros(Float64,L);
elastic = zeros(Float64,L);
drivingf = zeros(Float64,L);
aSites = zeros(Int64,L);
aSitesNew = zeros(Int64,L);
isActNew = zeros(Bool,L);
S_hist = zeros(0)
time_idx = 1:100000
for t in time_idx
    origin, is_relaxing, amount  = trigger!(thresholds, viscoel, elastic);
    S = avalanche!(origin, aSites, aSitesNew, isActNew, thresholds, viscoel, elastic, k0, k1, k2);
    append!(S_hist, S);
end
binning = 10 .^(range(0.0,5.0, length=50) )
histogram(S_hist, bins =binning,  norm=true, xaxis=(:log10, (1, 1e5)), yaxis=:log)
plot!(binning, 0.5*binning.^(-1.11));
plot!(binning, 0.5*binning.^(-1.45));
savefig("stuff.png");