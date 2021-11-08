from openmmtools import multistate
from simtk import unit
from simtk import unit

def cal(leg):
	reporter = multistate.MultiStateReporter(storage=f'{leg}.nc', open_mode='r', checkpoint_interval=1)
	analyzer = multistate.ReplicaExchangeAnalyzer(reporter=reporter)

	Deltaf_ij, dDeltaf_ij = analyzer.get_free_energy()
	DeltaH_ij, dDeltaH_ij = analyzer.get_enthalpy()

	# reference states in a simple REMD are 0 and -1 meaning the first and the last state are taken as the diff
	data = {}
	data['free_energy_diff'] = Deltaf_ij[0, -1]
	free_energy_diff_error = dDeltaf_ij[0, -1]

	phase_delta_f = data['free_energy_diff']
	# TODO: long range correction
	# phase_delta_f_ssc = data['free_energy_diff_standard_state_correction']

	phase_exp_delta_f = phase_delta_f # + phase_delta_f_ssc
	delta_f_unit = phase_exp_delta_f * analyzer.kT

	# per mole
	permole = delta_f_unit / unit.kilocalories_per_mole
	print(f'{leg} kcal/mol: {permole:.3f} ({free_energy_diff_error:.3f})')
	return permole, free_energy_diff_error


vac, verr = cal('vac')
sol, serr = cal('sol')
print(f'{vac - sol:.3f} ({(verr + serr)/2:.3f}), legs: vac {vac:.3f} ({verr:.3f}), sol {sol:.3f} ({serr:.3f}) kcal/mol')