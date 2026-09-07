import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from benchmark_conductive_cusp import run_benchmark
if __name__=='__main__':
 r=run_benchmark()
 assert all(r['gates'].values()),r
 print('PASS: manufactured conductive cusp transport and mesh refinement')
