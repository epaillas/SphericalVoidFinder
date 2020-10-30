import click
from python_tools.sphericalvoids import SphericalVoids


@click.command()
@click.option('--centres_filename', type=str, required=True)
@click.option('--tracers_filename', type=str, required=True)
@click.option('--ouput_filename', type=str, required=True)
@click.option('--ncores', type=int, default=1, help='Number of cores to use for parallel tasks.')
@click.option('--steps', type=str, default='1,2,3,4', help='Which steps are to be run. (e.g. 1,2,3).')
@click.option('--pos_cols', type=str, default='0,1,2', help='Indices of columns where tracer positions are stored.')
@click.option('--rvoidmax', type=float, default=100, help='Maximum void radius to search.')
@click.option('--box_size', type=float, default=1024, help='[Periodic box] Size of the simulation box')
@click.option('--ngrid', type=int, default=100)
def run_spherical_voids(centres_filename,
                        tracers_filename,
                        output_handle,
                        ncores,
                        steps,
                        pos_cols,
                        rvoidmax,
                        box_size,
                        ngrid):

    voids = SphericalVoids(centres_filename=centres_filename,
                           tracers_filename=tracers_filename,
                           output_handle=output_handle,
                           ncores=ncores,
                           steps=steps,
                           pos_cols=pos_cols,
                           rvoidmax=rvoidmax,
                           box_size=box_size,
                           ngrid=ngrid)


if __name__ == '__main__':
    run_spherical_voids()
