
# A gnuplot script that will plot relevant data from each subdirectory

# general settings for the run
set fit errorvariables

#vars
num_data = 8
#num_B30 = 8
#num_B100 = 7
array wrms[num_data]
array werr[num_data]
array wturb[num_data]
array wterr[num_data]
array wlam[num_data]
array wlerr[num_data]
array Urms[num_data]
array Uerr[num_data]
array urms[num_data]
array uerr[num_data]
array vrms[num_data]
array verr[num_data]
array tflux[num_data]
array tferr[num_data]
array Dmom[num_data]
array dmerr[num_data]
array Dtemp[num_data]
array dterr[num_data]
array power[num_data]
array perr[num_data]
array B[num_data]
array Pe[num_data]
array Re[num_data]
array invRo[num_data]
array files[num_data]
array lwb[num_data]
array upb[num_data]
array mix[num_data]
array mixerr[num_data]
array UHrms[num_data]
array UHerr[num_data]
array vortratio[num_data]
array vortratioerr[num_data]
array loop_per_file[num_data]
array dt[num_data]

# B30 simulations

files[1] = "<cat B30Re600Pe60/OUT*"
B[1] = 30.
invRo[1] = 0.
Re[1] = 600.
Pe[1] = 60.
lwb[1] = 38
upb[1] = 81


files[2] = "<cat Om0.5B30Re600Pe60/OUT*"
B[2] = 30.
invRo[2] = 0.5
Re[2] = 600.
Pe[2] = 60.
lwb[2] = 125.
upb[2] = 215.

files[3] = "<cat Om1B30Re600Pe60/OUT*"
B[3] = 30.
invRo[3] = 1.0
Re[3] = 600.
Pe[3] = 60.
lwb[3] = 216
upb[3] = 258

files[4] = "<cat Om2B30Re600Pe60/OUT*"
B[4] = 30.
invRo[4] = 2.
Re[4] = 600.
Pe[4] = 60.
lwb[4] = 247
upb[4] = 313

files[5] = "<cat Om3B30Re600Pe60/OUT*"
B[5] = 30.
invRo[5] = 3.
Re[5] = 600.
Pe[5] = 60.
lwb[5] = 325
upb[5] = 418

files[6] = "<cat Om4B30Re600Pe60/OUT*"
B[6] = 30.
invRo[6] = 4.
Re[6] = 600.
Pe[6] = 60.
lwb[6] = 429
upb[6] = 519

files[7] = "<cat Om5B30Re600Pe60/OUT*"
B[7] = 30.
invRo[7] = 5.
Re[7] = 600.
Pe[7] = 60.
lwb[7] = 485
upb[7] = 530

files[8] = "<cat Om10B30Re600Pe60/OUT*"
B[8] = 30.
invRo[8] = 10.
Re[8] = 600.
Pe[8] = 60.
lwb[8] = 534
upb[8] = 638

# Do loop to iterate over files and extract avg data
do for [i=1:num_data] {
    # U rms
    f(x) = a
    fit [lwb[i]:upb[i]] f(x) files[i] u 2:4 via a
    Urms[i] = a
    Uerr[i] = FIT_STDFIT
    loop_per_file[i] = int((upb[i]-lwb[i])*Urms[i]/20.)
    dt = 20./Urms[i]
    # u rms
    fit [lwb[i]:upb[i]] f(x) files[i] u 2:36 via a
    urms[i] = a
    uerr[i] = FIT_STDFIT
    # v rms
    fit [lwb[i]:upb[i]] f(x) files[i] u 2:37 via a
    vrms[i] = a
    verr[i] = FIT_STDFIT
    # w rms
    fit [lwb[i]:upb[i]] f(x) files[i] u 2:38 via a
    wrms[i] = a
    werr[i] = FIT_STDFIT
    # w_turb rms
    fit [lwb[i]:upb[i]] f(x) files[i] u 2:50 via a
    wturb[i] = a
    wterr[i] = FIT_STDFIT
    # w_lam rms
    fit [lwb[i]:upb[i]] f(x) files[i] u 2:51 via a
    wlam[i] = a
    wlerr[i] = FIT_STDFIT
    # t_flux
    fit [lwb[i]:upb[i]] f(x) files[i] u 2:9 via a
    tflux[i] = a
    tferr[i] = FIT_STDFIT
    # Dmom
    fit [lwb[i]:upb[i]] f(x) files[i] u 2:47 via a
    Dmom[i] = a
    dmerr[i] = FIT_STDFIT
    # Dtemp
    fit [lwb[i]:upb[i]] f(x) files[i] u 2:45 via a
    Dtemp[i] = a
    dterr[i] = FIT_STDFIT
    # Power 
    fit [lwb[i]:upb[i]] f(x) files[i] u 2:49 via a
    power[i] = a
    perr[i] = FIT_STDFIT
    fit [lwb[i]:upb[i]] f(x) files[i] u 2:((B[i]*$45/Pe[i])/(B[i]*$45/Pe[i]+$47/Re[i])) via a
    mix[i] = a
    mixerr[i] = FIT_STDFIT
    fit [lwb[i]:upb[i]] f(x) files[i] u 2:(($36**2 + $37**2)**(0.5)) via a
    UHrms[i] = a
    UHerr[i] = FIT_STDFIT
    fit [lwb[i]:upb[i]] f(x) files[i] u 2:(sqrt(($41**2)/($39**2 + $40**2 + $41**2))) via a
    vortratio[i] = a
    vortratioerr[i] = FIT_STDFIT
   
}

stats loop_per_file
num_loops = STATS_max

array Om0.5B30_avg_tflux[num_data, num_loops]
array Om0.5B30_avg_tflux_err[num_data, num_loops]
array Om0.5eff_invRo[num_data, num_loops]
array Om0.5eff_invRo_err[num_data, num_loops]
do for [i=1:num_data] {
    do for [j=1:loop_per_file[i]] {
        fit [lwb[i]:lwb[i]+j*dt] f(x) files[i] u 2:9 via a
        turn_avg_tflux[i][j] = a
        turn_avg_tflux_err[i][j] = FIT_STDFIT
        fit [lwb[i]:lwb[i]+j*dt] f(x) files[i] u 2:4 via a
        eff_invRo[i][j] = a*invRo
        eff_invRo_err[i][j] = FIT_STDFIT*invRo
    }
}


# another for loop to fin[i] u 2:9 wia a
#print B
save variables 'rot_eddy_avg.dat'
#print Urms
#print Uerr
#print wrms
#print werr
#print werr

