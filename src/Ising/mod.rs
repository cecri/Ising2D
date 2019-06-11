pub type Index = u32;

pub trait Model {
    fn neighbors(&self, n: Index) -> Vec<Index>;
    fn size(&self) -> u32;
    fn energy(&self, conf: &Vec<i8>) -> f64 {
        let mut res: f64 = 0.0;

        let all_neighbors = self.all_neighbors();
        for (i,j) in &all_neighbors {
            res -= f64::from(conf[*i as usize]*conf[*j as usize]);
        }
        res
    }
    fn magnetization(&self, conf: &Vec<i8>) -> f64;
    fn all_neighbors(&self) -> Vec<(Index,Index)>;
}


//n1 (row) X n2 (col), idx are counters as col major
pub struct Ising2D {
    pub n1: u32,
    pub n2: u32,
}

impl Ising2D {
    #[inline]
    pub fn to_idx(&self, i: u32, j: u32) -> u32 {
        j * self.n1 + i
    }
    #[inline]
    pub fn to_coord(&self, n: u32) -> (u32,u32) {
        (n%self.n1, n/self.n1)
    }
}

pub fn random_conf(n : u32) -> Vec<i8> {
    use rand::distributions::Uniform;
    use rand::prelude::*;

    let mut rng = rand::thread_rng();
    let dist = Uniform::from(0..2);
    let mut conf : Vec<i8> = Vec::with_capacity(n as usize);
    for _i in 0..n {
        conf.push( 1 - 2*dist.sample(&mut rng) );
    }
    conf
}

pub struct WolffSampler<'a, T: Model> {
    model: &'a T,
    beta: f64,
    conf: Vec<i8>,
}

impl Model for Ising2D{
    #[inline]
    fn size(&self) -> u32 {
        self.n1 * self.n2
    }

    fn neighbors(&self, n : Index) -> Vec<Index> {
        let (i,j) = self.to_coord(n);
        let mut res : Vec<u32> = Vec::new();
        if i >= 1 {
            res.push(self.to_idx(i-1,j));
        }
        if i+1 < self.n1 {
            res.push(self.to_idx(i+1,j));
        }
        if j >= 1 {
            res.push(self.to_idx(i,j-1));
        }
        if j +1 < self.n2 {
            res.push(self.to_idx(i,j+1));
        }

        res
    }

    fn magnetization(&self, conf: &Vec<i8>) -> f64 {
        let mut res: i32 = 0;
        for c in conf {
            res += *c as i32;
        }
        f64::from(res)
    }
    /*

    fn energy(&self, conf: &Vec<i8>) -> f64 {
        let mut res: f64 = 0.0;
        for i in 0..self.n1-1 {
            for j in 0..self.n2 {
                let n = i+j*self.n1;
                res -= f64::from(conf[n as usize]*conf[(n+1) as usize]);
            }
        }
        for i in 0..self.n1 {
            for j in 0..self.n2-1 {
                let n = i+j*self.n1;
                res -= f64::from(conf[n as usize]*conf[(n + self.n1) as usize]);
            }
        }
        res
    }*/

    fn all_neighbors(&self) -> Vec<(Index,Index)> {
        let mut res : Vec<(Index,Index)> = Vec::new();
        for i in 0..self.n1 {
            for j in 0..self.n2-1 {
                let n = j*self.n1 + i;
                res.push((n,n+self.n1));
            }
        }
        for i in 0..self.n1-1 {
            for j in 0..self.n2 {
                let n = j*self.n1 + i;
                res.push((n,n+1));
            }
        }
        res
    }
}

impl<'a, T: Model> WolffSampler<'a, T> {
    pub fn new (model: &'a T, beta: f64) -> WolffSampler <'a, T> {
        WolffSampler {
            model,
            beta,
            conf : random_conf(model.size()), 
        }
    }

    pub fn randomize_conf(&mut self) {
        self.conf = random_conf(self.size());
    }
    
    pub fn size(&self) -> u32 {
        self.model.size()
    }

    pub fn getConf(&self) -> &Vec<i8> {
        & self.conf
    }

    pub fn sweep (&mut self) {
        use std::collections::VecDeque;
        use std::collections::HashSet;
        use rand::distributions::Uniform;
        use rand::prelude::*;

        let mut cluster : HashSet<u32> = HashSet::new();
        let mut queue: VecDeque<u32> = VecDeque::new();


        let n = self.model.size();
        let idx_dist = Uniform::from(0..n);
        let p_dist : Uniform<f64> = Uniform::from(0.0..1.0);
        let mut rng = rand::thread_rng();

        let prob = 1.0 - (-2.0*self.beta).exp();

        let k = idx_dist.sample(&mut rng);
        queue.push_back(k);
        
        while let Some(k) = queue.pop_front() {
            cluster.insert(k);
            for x in self.model.neighbors(k) {
                if self.conf[k as usize] == self.conf[x as usize] {
                    //push x with probability p = 1-exp(-2*beta)
                    if !cluster.contains(&x) &&
                        (p_dist.sample(&mut rng) < prob) {
                            queue.push_back(x);
                    }
                }
            }
        }

        for x in cluster {
            self.conf[x as usize] *= -1;
        }
    }
}
