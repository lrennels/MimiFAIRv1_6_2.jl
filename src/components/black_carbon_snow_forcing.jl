# --------------------------------------------------
# Radiative forcing from black carbon on snow.
# --------------------------------------------------

@defcomp bc_snow_forcing begin

	E_ref_bc           = Parameter() # Reference year emissions.
	F_ref_bc           = Parameter() # Reference forcing value (W/m⁻²).
    BC_emiss        = Parameter(index=[time]) # Black carbon emissions (Mt yr⁻¹).

    forcing_BC_snow = Variable(index=[time])  # Radiative forcing from black carbon (Wm⁻²).


    function run_timestep(p, v, d, t)

   		# Caluclate forcing for BC.
	    # Assumes a scaling factor so the 2011 relationship between BC forcing (best AR5 estimate) and emissions (Meinshausen) holds for all years.
    	v.forcing_BC_snow[t] = p.BC_emiss[t] * (p.F_ref_bc / p.E_ref_bc)
	end
end
