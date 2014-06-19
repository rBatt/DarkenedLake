
irrInt <- function(I0, z1, z2, kd){
	# I0 = PAR at surface (above water)
	# z1 = shallow depth
	# z2 = deep depth
	# kd = light attenuation cofficient (slope of lm(light%~log(z)) when z ranges between z1 and z2)
	(-exp(-kd*z2)*I0)/kd - (-exp(-kd*z1)*I0)/kd # returns the average PAR between z1 and z2
}
irrInt(I0=1200, z1=0, z2=2, kd=1.15)



ifunc <- function(z, kd, i0) -exp(z*-kd)*i0

integrate(ifunc, upper=0, lower=2, kd=1.15, i0=1200)