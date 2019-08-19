extern crate nalgebra as na;
use super::Ising;

pub mod Utility
{
    use std::ops::{ShrAssign, BitAnd};
    pub fn to_sigma(n: u32, mut k: u32) -> Vec<i8> {
        let mut res : Vec<i8> = Vec::new();
        for i in 0..n {
            let t = (k & 1) as i8;
            res.push((1-2*t) as i8);
            k >>= 1;
        }
        res
    }
    pub fn to_int(conf: &Vec<i8>) -> u32 {
        let mut res : u32 = 0;
        for i in 0..conf.len() {
            res |= (1 << i)*(((1+conf[i])/2) as u32);
        }
        res
    }
}

pub trait Machine {
    fn get_num_visible(&self) -> u32;
    fn get_dim(&self) -> u32;
    fn coeff(&self, k : u32) -> f64;
    fn partial_der(&self, sigma : &na::DVector<i8>)  -> na::DVector<f64>;

    fn from_ising<T: Ising::Model>(ising: &T, beta: f64) -> Self;
}

pub struct RBM {
    n: u32,
    m: u32,
    pub weights: na::DMatrix<f64>,
}

unsafe impl Send for RBM {}
unsafe impl Sync for RBM {}


impl Machine for RBM {
    fn get_num_visible(&self) -> u32 {
        self.n
    }



    fn get_dim(&self) -> u32 {
        self.n * self.m
    }

    fn coeff(&self, k: u32) -> f64 {
        use na::{DMatrix, DVector};
        use Utility::to_sigma;

        let s : DVector<f64> = na::convert(DVector::from_vec(to_sigma(self.n,k)));

        let theta = &self.weights*&s;
        let mut res:f64 = 1.0;
        for t in &theta {
            res *= 2.0*t.cosh();
        }
        res
    }

    fn from_ising<T: Ising::Model>(ising: &T, beta: f64) -> RBM
    {
        use na::{DMatrix, DVector};

        let n = ising.size();
        let neighbors = ising.all_neighbors();
        let m = neighbors.len();
        let mut weights: DMatrix<f64> = DMatrix::zeros(m, n as usize);
        let w = (beta).exp().acosh()/2.0;

        for (count, pair) in neighbors.iter().enumerate() {
            weights[(count,pair.0 as usize)] = w;
            weights[(count,pair.1 as usize)] = w;
        }
        RBM{
            n,
            m: m as u32,
            weights
        }
    }
    
    fn partial_der(&self, sigma : &na::DVector<i8>)  -> na::DVector<f64> {
        use na::DVector;
        let s : DVector<f64> = na::convert_ref(sigma);
        let mut res: DVector<f64> = DVector::zeros((self.n*self.m) as usize);
        let theta = &self.weights*&s;
        for i in 0..self.n {
            for j in 0..self.m {
                res[(i*self.m + j) as usize] = f64::from(sigma[i as usize])*theta[j as usize].tanh();
            }
        }
        res
    }


}

impl RBM {
    pub fn get_num_hidden(&self) -> u32 {
        self.m
    }
    pub fn get_dim_full(&self) -> u32 {
        self.n + self.m + self.n*self.m
    }
    pub fn partial_der_full(&self, sigma : &na::DVector<i8>)  -> na::DVector<f64> {
        use na::DVector;
        let s : DVector<f64> = na::convert_ref(sigma);
        let mut res: DVector<f64> = DVector::zeros((self.n + self.m + self.n*self.m) as usize);
        let theta = &self.weights*&s;
        for i in 0..self.n {
            res[i as usize] = f64::from(sigma[i as usize]);
        }
        for j in 0..self.m {
            res[(self.n + j) as usize] = theta[j as usize].tanh();
        }

        for i in 0..self.n {
            for j in 0..self.m {
                let idx = i*self.m + j;
                res[(self.n + self.m + idx) as usize] = f64::from(sigma[i as usize])*theta[j as usize].tanh();
            }
        }
        res
    }
}


pub fn coeffs<T: Machine>(machine: &T, normalize: bool) -> na::DVector<f64> {
    let dim = 1 << machine.get_num_visible();
    let mut c : na::DVector<f64> = na::DVector::zeros(dim);
    for i in 0..dim {
        c[i] = machine.coeff(i as u32);
    }
    if normalize {
        c.normalize_mut();
    }
    c
}
