"""
Manages a pool of processtex jobs.

Don't do this in processtex, to avoid crashes...
"""

import argparse
import glob
import os

from multiprocessing import Pool, cpu_count
from random import shuffle
from subprocess import Popen


PROCESSTEX = os.path.join(os.path.dirname(__file__), 'processtex.py')


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i+n]


def job(arg):
    """Run one processing job."""
    args, htmls = arg
    cmdline = [
        'python3', PROCESSTEX,
        '--preamble', args.preamble,
        '--style-path', args.style_path,
        '--cache-dir', args.cache_dir,
        '--img-dir', args.img_dir,
        '--output-dir', args.output_dir,
        '--build-dir', args.build_dir,
    ]
    if args.no_cache:
        cmdline.append('--no-cache')
    proc = Popen(cmdline + htmls)
    proc.wait()
    if proc.returncode != 0:
        raise Exception("Call failed")


def main():
    """Run the main function."""
    parser = argparse.ArgumentParser(
        description='Process LaTeX in html files: job dispatcher.')
    parser.add_argument('--preamble', default='preamble.tex', type=str,
                        help='LaTeX preamble')
    parser.add_argument('--style-path', default='', type=str,
                        help='Location of LaTeX style files')
    parser.add_argument('--cache-dir', default='pretex-cache', type=str,
                        help='Cache directory')
    parser.add_argument('--img-dir', default='figure-images', type=str,
                        help='LaTeX image include directory')
    parser.add_argument('--no-cache', action='store_true',
                        help='Ignore cache and regenerate')
    parser.add_argument('--chunk-size', type=int, default=50,
                        help='Run processtex on chunks of this size')
    parser.add_argument('--output-dir', type=str, required=True,
                        help='HTML output directory')
    parser.add_argument('--build-dir', type=str, required=True,
                        help='Final build directory')
    args = parser.parse_args()

    os.makedirs(os.path.join(args.build_dir, 'knowl'), exist_ok=True)

    htmls = \
        glob.glob(os.path.join(args.output_dir, '*.html')) + \
        glob.glob(os.path.join(args.output_dir, 'knowl', '*.html'))

    # Process in a random order.  Otherwise one process gets all the section
    # files.
    shuffle(htmls)
    # htmls = [os.path.join(args.output_dir, 'similarity.html'),
    #          os.path.join(args.output_dir, 'index.html'),
    #          # os.path.join(args.output_dir, 'index2.html'),
    #          ]

    with Pool(processes=max(cpu_count()-1, 3)) as pool:
        job_args = []
        for chunk in chunks(htmls, args.chunk_size):
            job_args.append((args, chunk))
        result = pool.map_async(
            job, job_args, error_callback=lambda x: pool.close())
        result.wait()


if __name__ == "__main__":
    main()
