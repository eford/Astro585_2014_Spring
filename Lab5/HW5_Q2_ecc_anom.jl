# Calculate accentric annomaly from mean anommaly and eccentricity [to within tol]
# Based on algorithm from Danby
function ecc_anom(mean_anom::Real, ecc::Real; tol::Real = 1.e-12)
    #M = mod2pi(mean_anom)
    M = mod(mean_anom,2pi)
    if(M<0.) M += 2pi end
    @assert(0. <= M <= 2pi) # should be redundant
    @assert(0. <= ecc < 1.0)
    @assert(0 < tol <= 1.0e-4)
    const num_max_it = 20
    const k = 0.85
    const third = 1.0/3.0
    x = (M<pi) ? M + k*ecc : M - k*ecc;
	es = ecc*sin(x)
    ec = ecc*cos(x)
    F = (x-es)-M;
    count = 1
    while abs(F)>tol 
	  Fp = 1.-ec;
	  Fpp = es;
	  Fppp = ec;
	  Dx = -F/Fp;
	  Dx = -F/(Fp+0.5*Dx*Fpp);
	  Dx = -F/(Fp+0.5*Dx*(Fpp+third*Dx*Fppp));
	  x += Dx;
	  es = ecc*sin(x)
      ec = ecc*cos(x)
      F = (x-es)-M;
      count += 1
      if count>num_max_it error("ecc_anom($mean_anom, $e, $tol) did not converge!\n") end
    end
  x
end

# Calculate eccentric annomaly for array of mean anommalies and one eccentricity [to within tol]
function ecc_anom(mean_anom::Array{Float64}, ecc::Float64; tol::Real = 1.e-10)
  E = similar(mean_anom)
  for i in 1:length(mean_anom)
      E[i] = ecc_anom(mean_anom[i],ecc,tol=tol)
  end
  E      
end

# Calculate accentric annomaly for array of mean anommalies and an array of eccentricities [to within tol]
function ecc_anom(mean_anom::Array{Float64}, ecc::Array{Float64}; tol::Real = 1.e-10)
  @assert(size(mean_anom)==size(ecc))
  E = similar(mean_anom)
  for i in 1:length(mean_anom)
      E[i] = ecc_anom(mean_anom[i],ecc[i],tol=tol)
  end
  E      
end

function decc_anom_dM(mean_anom::Real, ecc::Real; tol::Real = 1.e-10)
  1.0/(1.0-ecc*ecc_anom(mean_anom,ecc,tol=tol))
end
