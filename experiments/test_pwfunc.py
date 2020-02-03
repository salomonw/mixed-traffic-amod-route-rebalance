import src.CARS as cars

fcoeffs = [1,0,0,0,0.15,0]

fun = cars.get_approx_fun(fcoeffs, nlines=2)

fun = cars.get_approx_fun(fcoeffs, xa=1.2, nlines=2)