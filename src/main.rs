use std::vec;


mod Ising2D {
    pub struct Model {
        pub n1: u32,
        pub n2: u32,
    }
    pub struct WolffSampler<'a> {
        model: &'a Model,
        beta: f64,
        conf: Vec<i8>,
    }

    impl Model {
        pub fn size(&self) -> u32 {
            self.n1 * self.n2
        }
    }

    impl<'a> WolffSampler<'a> {
        pub fn new (model: &'a Model, beta: f64) -> WolffSampler <'a> {
            use rand::distributions::Uniform;
            use rand::prelude::*;

            let dist = Uniform::from(0..2);
            let mut rng = rand::thread_rng();

            let n = (model.n1*model.n2) as usize;
            let mut conf : Vec<i8> = Vec::with_capacity(n);

            for _i in 0..n {
                conf.push( 1 - 2*dist.sample(&mut rng) );
            }

            WolffSampler {
                model,
                beta,
                conf, 
            }
        }
        
        pub fn size(&self) -> u32 {
            self.model.size()
        }

        pub fn getConf(&self) -> &Vec<i8> {
            & self.conf
        }
        
        fn neighbors(&self, n : u32) -> Vec<u32> {
            let mut res : Vec<u32> = Vec::new();
            res.push(n + self.model.n1);
            if n > self.model.n1 {
                res.push(n - self.model.n1);
            }
            res.push(n + 1);
            if n > 1 {
                res.push(n - 1);
            }
            
            res.retain(|&i| i < self.size());
            res
        }

        pub fn sweep (&mut self) {
            use std::collections::VecDeque;
            use rand::distributions::Uniform;
            use rand::prelude::*;

            let mut cluster : Vec<u32> = Vec::new();
            let mut queue: VecDeque<u32> = VecDeque::new();


            let n = self.model.n1*self.model.n2;
            let idx_dist = Uniform::from(0..n);
            let p_dist : Uniform<f64> = Uniform::from(0.0..1.0);
            let mut rng = rand::thread_rng();

            let prob = 1.0 - self.beta.exp();
            
            let k = idx_dist.sample(&mut rng);
            queue.push_back(k);
            println!("{}", k);
            
            while let Some(k) = queue.pop_front() {
                cluster.push(k);
                for x in self.neighbors(k) {
                    if self.conf[k as usize] == self.conf[x as usize] {
                        //push x with probability p = 1-exp(-2*beta)
                        if p_dist.sample(&mut rng) < prob
                        {
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
}

fn main() {
    use Ising2D;
    let ising2D = Ising2D::Model {
        n1: 5,
        n2: 5,
    };

    let mut ws = Ising2D::WolffSampler::new(&ising2D, 0.1);

    println!("{:?}", ws.getConf());
    ws.sweep();
    println!("{:?}", ws.getConf());
}
