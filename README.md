## Example cylinder telescope simulation

This directory contains the ingredients to run a sample simulation on Cedar, including generating telescope products (beam transfer matrices, telescope-SVD basis vectors, double-KL basis vectors, and power spectrum Fisher matrices) and simulated data. The telescope configuration is Simon's "half-pathfinder" (e.g. [doclib #1256](https://bao.chimenet.ca/doc/documents/1256)).

To proceed, copy this directory to your scratch or project space on Cedar, and make sure you can load a CHIME environment with the standard [radiocosmology](https://github.com/radiocosmology/) packages (`caput`, `cora`, `draco`, `driftscan`). You'll have to edit file paths in some files, as indicated below. (This could be done more intelligently using environment variables...but I'll leave that as an exercise!)

Please let Simon know if you encounter any issues!

### Generating telescope products

1. In `bt.yaml`, edit `output_directory` so that it points to your own `ch_scripts/sim_example/sim_output/products/` directory, and also edit `venv` to point to your Python virtual environment directory. From the command line, execute `slurm_generate_telescope_products.sh`. This uses the `drift-makeproducts` script in `driftscan` to submit a job to cedar's batch queue that generates telescope beam transfer matrices, SVD compressions of those matrices, KL basis vectors, and power spectrum Fisher matrices, using the config in `bt.yaml`.

### Generating simulated data

1. Submit an interactive job to cedar by running, e.g., `salloc --time=1:00:0 --nodes=1 --ntasks-per-node=1 --cpus-per-task=48 --mem=0 --account=rpp-chime --job-name=maps`. This will log you into a cedar compute node. From here, activate the standard CHIME modules and your own virtual environment, then run `generate_input_maps.sh` from the command line. This will use `cora` to generate simulated 21cm and foreground maps inside `sim_output/cora_maps/`. This should only take a few minutes, and you can exit the interactive job afterwards by logging out of your node.

2. In `sstream.yaml`, edit the venv location and the path to your simulation output directory, then run `slurm_gen_sstreams.sh` from the command line. This uses the `caput-pipeline` script in `caput` to submit a job to cedar's batch queue that computes simulated sidereal streams based on the sky maps generated above, using the config in `sstream.yaml`.

3. In `meas.yaml`, edit the venv location and the path to your simulation output directory, then run `slurm_measure.sh` from the command line. This submits a batch job that measures m-modes, telescope-SVD modes, KL modes, and power spectrum bandpowers from the simulated sidereal stream with noise, using the config in `meas.yaml`.

The notebook `plotting_power_spectrum.ipynb` demonstrates how to use routines in `ps_plot_utils.py` to visualize the Fisher errorbars on the power spectrum, along with the power spectrum measured from the simulated data.


### Generating simulated data corrupted by beam width perturbations

1. In `bt_pert.yaml`, edit `output_directory` so that it points to your own `ch_scripts/sim_example/sim_output/products_pert/` directory, and also edit `venv` to point to your Python virtual environment directory. From the command line, execute `slurm_generate_perturbed_telescope_products.sh`.

2. In `slurm_gen_sstream_pert_step1.sh`, edit the venv location and the path to your simulation output directory, then run `slurm_gen_sstream_pert_step1.sh` from the command line. This submits a batch job that simulates a sidereal stream that separately includes an unperturbed piece and a piece corresponding to beam-width perturbations, using the config in `sstream_pert_step1.yaml`.

3. Make the same edits to `slurm_gen_sstream_pert_step2.sh` and run it from the command line. This submits a batch job that generates a sidereal stream in the presence of beam-width perturbations with rms 0.1%, using the config in `sstream_pert_step2_std0.001.yaml`. IMPORTANT: The routine for this is not yet merged into the `draco` master branch (as of October 2021), so you need to checkout the `sf/pert-beam` branch for this step (and only this step).

4. Make the same edits to `slurm_measure_pert.sh` and run it from the command line. This makes the same measurements as `slurm_measure.sh` but using the perturbed sidereal stream, using the config in `meas_pert.yaml`.


## No-SVD simulation

Everything above uses the combination of the telescope-SVD and KL filtering in the power spectrum measurement process, but we can also skip the SVD step and rely solely on the KL filter. Currently (as of 20 October 2021), the no-SVD option is broken in the master branch of `driftscan`, so you'll need to checkout the `sf/nosvd-fix` branch to run these steps.

1. Edit the file paths in `make_beam_m_symlink.sh` and then run it from the command line. This will create new directories to hold the no-SVD telescope products, and also create a symlink to the previously-computed beam transfer matrices (stored by default in `sim_output/products/bt/beam_m`). The beam transfer matrices are identical whether or not the telescope-SVD is used, so the symlink prevents us from having to needlessly recompute them.

2. Edit the `venv` and `output_directory` in `bt_nosvd.yaml` consistently with the no-SVD output paths specified in `make_beam_m_symlink.sh`, and then run `slurm_generate_nosvd_telescope_products.sh`. This will run on the batch queue, and take a long time (around 8 hours).

3. Make appropriate edits to `meas.yaml`, and then run `slurm_measure_nosvd.sh`. This will measure the power spectrum of the previously-generated sidereal stream using the no-SVD KL filter and power spectrum pipeline. At this point, you can run the cells near the end of the `plotting_power_spectrum.ipynb` notebook to compare the output with and without the SVD step.

