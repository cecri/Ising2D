extern crate nalgebra as na;
extern crate itertools;
extern crate itertools_num;
extern crate nalgebra_lapack as nala;
extern crate scoped_threadpool;
extern crate serde;
extern crate serde_pickle;

mod Ising;
mod Machines;
use std::error::Error;
use Machines::{Machine, Utility, coeffs,RBM};

fn machine_exact_energy<T: Ising::Model, U: Machines::Machine>(ising: &T, rbm: &U) -> f64 {
    let cfs = coeffs(rbm, true);
    let mut energy: f64 = 0.0;

    for i in 0..cfs.len() {
        let conf = Utility::to_sigma(rbm.get_num_visible(), i as u32);
        let e = ising.energy(&conf);
        energy += cfs[i].powi(2)*e;
    }
    energy
}
fn ising_exact_prob<T: Ising::Model>(ising: &T, beta: f64) -> Vec<f64> {
    let n = ising.size();
    let dim = 1<<n;
    
    let mut prob: Vec<f64> = Vec::new();
    for i in 0..dim {
        let conf = Utility::to_sigma(n, i as u32);
        let e = ising.energy(&conf);
        prob.push((-beta*e).exp());
    }
    let s: f64 = prob.iter().sum();

    let res: Vec<f64> = prob.iter().map(|x| x/s).collect();
    res
}

fn sampled_energy<T: Ising::Model>(model: &T, samples: &Vec<Vec<i8>>) -> f64 {
    let egs : Vec<f64> = samples.iter().map(|x| model.energy(x) ).collect();
    egs.iter().sum::<f64>() / (egs.len() as f64)
}
fn sampled_mag<T: Ising::Model>(model: &T, samples: &Vec<Vec<i8>>) -> f64 {
    let mgs : Vec<f64> = samples.iter().map(|x| model.magnetization(x).abs() ).collect();
    mgs.iter().sum::<f64>() / (mgs.len() as f64)
}

fn ising_sample<T: Ising::Model>(model: &T, beta: f64, nsmp: usize) -> Vec<Vec<i8> > {
    let mut res: Vec<Vec<i8> > = Vec::new();
    let mut ws = Ising::WolffSampler::new(model, beta);

    for _j in 0..100 {
        ws.sweep();
    }

    for _i in 0..nsmp {
        res.push(ws.getConf().to_vec());
        ws.sweep();
    }
    res
}

struct SmatrixBuilder {
    pub g2 : na::DMatrix<f64>,
    pub g : na::DVector<f64>,
}

fn corr_fn(idx0: usize, idx1: usize, samples: &Vec<Vec<i8>>) -> f64
{
    let nsmp = samples.len() as f64;
    let two_points: i32 = samples.iter().map(|x| (x[idx0]*x[idx1]) as i32).sum();
    let one_points0: i32 = samples.iter().map(|x| x[idx0] as i32).sum();
    let one_points1: i32 = samples.iter().map(|x| x[idx1] as i32).sum();

    f64::from(two_points)/nsmp - f64::from(one_points0)/nsmp*f64::from(one_points1)/nsmp
}

impl SmatrixBuilder {
    fn construct_from_samples<T: Machine + Sync>(machine: &T, samples: &Vec<Vec<i8> >) -> SmatrixBuilder {
        use na::{DMatrix, DVector};
        use std::ops::AddAssign;
        use scoped_threadpool::Pool;
        use std::ptr::copy_nonoverlapping;
        
        let dim = machine.get_dim() as usize;
        let n_smp = samples.len();
        let mut gs: DMatrix<f64> = DMatrix::zeros(dim, n_smp);
        
        let mut pool = Pool::new(4);
        
        pool.scoped(|scoped| {
            for (i, mut col) in gs.column_iter_mut().enumerate(){
                scoped.execute(move || {
                    let sigma = na::DVector::<i8>::from_vec(samples[i].clone());
                    let src = machine.partial_der(&sigma);
                    col += src;
                });
            } //for
        });

        let mut g2 = &gs*gs.transpose();
        g2 /= f64::from(n_smp as u32);

        SmatrixBuilder {
            g2,
            g : gs.column_mean()
        }
    }
    fn construct_full_from_samples(machine: &RBM, samples: &Vec<Vec<i8> >) -> SmatrixBuilder {
        use na::{DMatrix, DVector};
        use std::ops::AddAssign;
        use scoped_threadpool::Pool;
        use std::ptr::copy_nonoverlapping;
        
        let dim = machine.get_dim_full() as usize;
        let n_smp = samples.len();
        let mut gs: DMatrix<f64> = DMatrix::zeros(dim, n_smp);
        
        let mut pool = Pool::new(4);
        
        pool.scoped(|scoped| {
            for (i, mut col) in gs.column_iter_mut().enumerate(){
                scoped.execute(move || {
                    let sigma = na::DVector::<i8>::from_vec(samples[i].clone());
                    let src = machine.partial_der_full(&sigma);
                    col += src;
                });
            } //for
        });

        let mut g2 = &gs*gs.transpose();
        g2 /= f64::from(n_smp as u32);

        SmatrixBuilder {
            g2,
            g : gs.column_mean()
        }
    }

    fn get_smat(&self) -> na::DMatrix<f64> {
        let mut res = self.g2.clone();
        res -= &self.g*self.g.transpose();
        res
    }
}

fn ent_ee(n: usize, m: usize, vec: &na::DVector<f64>) -> f64 {
    assert!(vec.len() == n*m);
    let mut psi = na::DMatrix::<f64>::zeros(m,n);
    psi.copy_from_slice(vec.as_slice());
    println!("{}\n", psi);

    let rho = &psi*psi.transpose();

    println!(r"Tr[\rho] = {}\n", rho.trace());
    let lambdas = nala::SymmetricEigen::eigenvalues(rho);
    
    let mut res: f64 = 0.0;
    for lambda in lambdas.into_iter() {
        if *lambda < 1e-10 {
            continue;
        }
        res -= lambda*lambda.log2();
    }
    res
}

fn main() -> Result<(), Box<dyn Error>>{
    use Ising::{Model, Ising2D};
    use Machines::{Machine, RBM, coeffs};
    use core::cmp::max;
    use std::fs::File;
    use std::io::Write;
    use na::DVector;

    
    let ising2D = Ising2D {
        n1: 8,
        n2: 8,
    };


    //let betas = vec![0.01, 0.30, 0.44, 0.58, 0.65, 0.9];
    //let betas = vec![2.0];
    let nsmp: usize = 10000;

    let betas = itertools_num::linspace(0.1, 1.0, 91);

    for beta in betas {
        let smps = ising_sample(&ising2D, beta, nsmp);

        println!("beta: {:.3}", beta);

        let rbm = RBM::from_ising(&ising2D, beta);
        let n = rbm.get_num_visible() as usize;
        let m = rbm.get_num_hidden() as usize;
        
        eprintln!("Constructing S Matrix Start");
        let smb = SmatrixBuilder::construct_full_from_samples(&rbm, &smps);
        let smat =  smb.get_smat();
        eprintln!("Constructing S Matrix Finished");
        eprintln!("Calculating eigenvalues Start");
        let evs = nala::SymmetricEigen::eigenvalues(smat);
        eprintln!("Calculating eigenvalues Finished");

        let filename = format!("EIGS_{:03}.dat", (beta*100.0 + 0.5) as u32);
        let mut file = File::create(filename).unwrap();

        //let length = eigs.eigenvalues.len();

        write!(&mut file, "{}\n", beta)?;
        for v in evs.iter() {
            write!(&mut file, "{}\n", v)?;
        }
    }

    Ok(())
}
