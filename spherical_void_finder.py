import click
from python_tools.sphericalvoids import SphericalVoids


@click.command()
@click.option('--centres_filename', type=str, required=True)
@click.option('--tracers_filename', type=str, required=True)
@click.option('--output_handle', type=str, required=True)
@click.option('--ncores', type=int, default=1)
@click.option('--steps', type=str, default='1,2,3,4')
@click.option('--pos_cols', type=str, default='0,1,2')
@click.option('--rvoidmax', type=float, default=100)
@click.option('--box_size', type=float, default=1024)
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
