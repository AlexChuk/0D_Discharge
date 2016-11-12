# include "0D-DischMod.h"
# include "Initials.h"
# include "EEDF_calc.h"
# include "Chemistry_calc.h"
# include "Gas_calc.h"

int main(void)
{
	init_read();
	init_gasDBparser();
	init_data();

	Nedf = EEDF_read_CS();
	Nchem = chem_make_react(Nedf);
	chem_read_react(Nchem);

	int nt,dot;
	Nt = int(tau/dt);
	dot = 0;
	for(nt=0;nt<=Nt;nt++)
	{
		dot += 1;

        EEDF_calc(nt,dot);
        EEDF_const_calc();

		chem_const(Te,Tgas,nt,dot);
		chem_runge_kutta4(Ni,nt,dot);

		gas_TP_calc(Ni,Tgas,Pgas,Hgas,nt,dot);

		if(dot==Ndots)
			dot = 0;
	}

	return 0;
}
