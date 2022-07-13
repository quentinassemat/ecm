use crate::primes::*;
use gmp::mpz::*;
use gmp::rand::*;
use std::fmt;
use std::fs::File;
use std::io::Write;
use std::time::{SystemTime, UNIX_EPOCH};

//Montgomery's elliptic curve arithmetic :

#[derive(Clone, Debug)]
pub struct Point {
    pub x: Mpz,
    pub z: Mpz,
}

impl Point {
    pub fn zero() -> Point {
        Point {
            x: Mpz::one(),
            z: Mpz::zero(),
        }
    }
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({} : {})", self.x, self.z)
    }
}

#[derive(Clone, Debug)]
pub struct Factor {
    pub prime: Mpz,
    pub exp: u32,
}

impl Factor {
    pub fn new(prime: Mpz) -> Factor {
        Factor { prime, exp: 0 }
    }

    pub fn incr(&mut self) {
        self.exp += 1;
    }
}

impl fmt::Display for Factor {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}^{}", self.prime, self.exp)
    }
}

pub fn print_factors(factors: &[Factor], f: &mut File) {
    if !factors.is_empty() {
        for item in factors.iter().take(factors.len() - 1) {
            let buf = format!("{} * ", item);
            f.write(buf.as_bytes()).expect("Erreur write");
            // print!("{} * ", item );
        }
        let buf = format!("{}\n", factors[factors.len() - 1]);
        f.write(buf.as_bytes()).expect("Erreur write");
        // println!("{}", factors[factors.len() - 1])
    } else {
        f.write(b"empty").expect("Erreur write");
        // println!("empty")
    }
}

pub fn x_add(p: &Point, q: &Point, pmq: &Point, n: &Mpz) -> Point {
    let u: Mpz = (&p.x - &p.z % n) * (&q.x + &q.z % n) % n;
    let v: Mpz = (&p.x + &p.z % n) * (&q.x - &q.z % n) % n;
    let upv = (&u + &v).powm(&Mpz::from(2), n);
    let umv = (&u - &v).powm(&Mpz::from(2), n);
    Point {
        x: &pmq.z * &upv % n,
        z: &pmq.x * &umv % n,
    }
}

pub fn x_dbl(p: &Point, n: &Mpz, a: &Mpz) -> Point {
    let q = (&p.x + &p.z).powm(&Mpz::from(2), n);
    let r = (&p.x - &p.z).powm(&Mpz::from(2), n);
    let s = &q - &r % n;
    Point {
        x: &q * &r % n,
        z: (&s * (&r + a * &s)) % n,
    }
}

// on a enlevé le constant
pub fn ladder(m: &Mpz, x: &Point, n: &Mpz, a: &Mpz) -> Point {
    let beta = m.size_in_base(2);
    let mut bit: u8;
    let mut x1 = x.clone();
    let mut x0 = Point::zero();
    for i in 0..beta {
        bit = m.tstbit(beta - i - 1) as u8; // bit vaut 1 ou 0
        if bit == 0 {
            let temp = x_dbl(&x0, n, a);
            x1 = x_add(&x0, &x1, x, n);
            x0 = temp;
        } else {
            let temp = x_dbl(&x1, n, a);
            x0 = x_add(&x0, &x1, x, n); // valide car x_add(&x0, &x1, &x) = x_add(&x1, &x0, &x)
            x1 = temp;
        }
    }
    x0
}

pub fn trial_division(n: &Mpz) -> Vec<Factor> {
    let mut res: Vec<Factor> = Vec::new();
    let mut prime = Mpz::from(2);
    let mut div_test = prime.clone();
    // on fait le cas 2 à part :
    {
        let mut f: Factor = Factor::new(prime.clone());
        while n.is_multiple_of(&div_test) {
            div_test *= &prime;
            f.incr();
        }
        if f.exp > 0 {
            res.push(f);
        }
        prime += 1;
        div_test.set(&prime);
    }
    for add in PRIMES {
        let mut f: Factor = Factor::new(prime.clone());
        while n.is_multiple_of(&div_test) {
            div_test *= &prime;
            f.incr();
        }
        if f.exp > 0 {
            res.push(f);
        }
        prime += (2 * add) as u64;
        div_test.set(&prime);
    }
    res
}

// sigma \in Z / nZ \{1, 2, 3, 4}
pub fn random_curve(sigma: &Mpz, n: &Mpz) -> (Mpz, Point) {
    let (u, v) = ((sigma * sigma - 5) % n, (4_u64 * sigma) % n);
    let three = Mpz::from(3);
    let p = Point {
        x: Mpz::powm(&u, &three, n),
        z: Mpz::powm(&v, &three, n),
    };
    let a = (Mpz::powm(&(&v - &u), &three, n)
        * (3_u64 * &u + &v)
        * Mpz::invert(&(&(4_u64 * Mpz::powm(&u, &three, n) * v) % n), n)
            .expect("Error inversion\n")
        - 2)
        % n;
    (a, p)
}

pub fn random_mpz(n: &Mpz) -> Mpz {
    let mut rng = RandState::new();
    let time: u64 = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos() as u64;
    rng.seed_ui(time);
    rng.urandom(&(n - 5)) + 5
}

// ECMTrial with only stage 1
#[warn(dead_code)]
pub fn ecm_trial(n: &Mpz, bound: &Mpz, f: &mut File) -> Option<Mpz> {
    let sigma = random_mpz(n);
    let (a, p) = random_curve(&sigma, n);
    let test = (&a * &a - 4).gcd(n);
    if test != Mpz::one() {
        let buf = format!("factor : {}, a : {}, p : {}\n", test, a, p);
        f.write(buf.as_bytes()).expect("Erreur write");
        // println!("factor : {}, a : {}, p : {}", test, a, p);
        return Some(test);
    }
    let mut q = p.clone();
    let mut l = Mpz::from(2);
    while l < *bound {
        let mut m = l.clone();
        while m < *bound {
            m *= &l;
        }
        m = m.div_floor(&l);
        q = ladder(&m, &q, n, &((&a + 2) / 4));
        l = l.nextprime();
    }
    let g = q.z.gcd(n);
    if Mpz::one() < g && g < *n {
        let buf = format!("factor : {}, a : {}, p : {}\n", g, a, p);
        f.write(buf.as_bytes()).expect("Erreur write");
        // println!("factor : {}, a : {}, p : {}", g, a, p);
        return Some(g);
    }
    None
}

// ECMTrial with stage 2
pub fn ecm_trial2(n: &Mpz, bound1: &Mpz, bound2: &Mpz, f: &mut File) -> Option<Mpz> {
    let sigma = random_mpz(n);
    let (a, p) = random_curve(&sigma, n);
    let test = (&a * &a - 4).gcd(n);
    if test != Mpz::one() {
        let buf = format!("factor : {}, a : {}, p : {}\n", test, a, p);
        f.write(buf.as_bytes()).expect("Erreur write");
        // println!("factor : {}, a : {}, p : {}", test, a, p);
        return Some(test);
    }
    let mut q = p.clone();
    let mut l = Mpz::from(2);
    while l < *bound1 {
        let mut m = l.clone();
        while m < *bound1 {
            m *= &l;
        }
        m = m.div_floor(&l);
        q = ladder(&m, &q, n, &((&a + 2) / 4));
        l = l.nextprime();
    }
    let mut g = q.z.gcd(n);
    if Mpz::one() < g && g < *n {
        let buf = format!("factor : {}, a : {}, p : {}\n", g, a, p);
        f.write(buf.as_bytes()).expect("Erreur write");
        // println!("factor : {}, a : {}, p : {}", g, a, p);
        return Some(g);
    }
    //stage 2 :
    while l < *bound2 {
        let q2 = ladder(&l, &q, n, &((&a + 2) / 4));
        g = q2.z.gcd(n);
        if Mpz::one() < g && g < *n {
            let buf = format!("factor : {}, a : {}, p : {}\n", g, a, p);
            f.write(buf.as_bytes()).expect("Erreur write");
            // println!("factor : {}, a : {}, p : {}", g, a, p);
            return Some(g);
        }
        l = l.nextprime();
    }
    None
}

pub fn factorization(n: &Mpz, bound1: &Mpz, bound2: &Mpz, f: &mut File) {
    let mut factors = trial_division(n);
    let mut work = n.clone();
    let buf = format!("Factors with trial_division\n");
    f.write(buf.as_bytes()).expect("Erreur write");
    // println!("Factors with trial_division");
    print_factors(&factors, f);
    // on enlève les petits facteurs de n
    for f in &factors {
        work = work.div_floor(&Mpz::powm(&f.prime, &Mpz::from(f.exp), n));
    }
    let buf = format!("Integer without trivial factors : {}\n", &work);
    f.write(buf.as_bytes()).expect("Erreur write");
    // println!("Integer without trivial factors : {}", &work);
    // on lance la recherche de facteur non trivial maintenant
    let buf = format!(
        "Lauching ECM with bound1 : {}, bound2 : {}\n",
        bound1, bound2
    );
    f.write(buf.as_bytes()).expect("Erreur write");
    // println!("Lauching ECM with bound1 : {}, bound2 : {}", bound1, bound2);
    while work.probab_prime(50) == ProbabPrimeResult::NotPrime {
        if let Some(f) = ecm_trial2(&work, bound1, bound2, f) {
            let mut div_test = f.clone();
            let mut fact: Factor = Factor::new(f.clone());
            while work.is_multiple_of(&div_test) {
                div_test *= &f;
                fact.incr();
            }
            work = work.div_floor(&Mpz::powm(&fact.prime, &Mpz::from(fact.exp), n));
            factors.push(fact);
        }
    }
    factors.push(Factor {
        prime: work,
        exp: 1,
    });
    f.write(b"n = ").expect("Erreur write");
    // print!("n = ");
    print_factors(&factors, f);
    f.write(b"End ECM\n").expect("Erreur write");
    // println!("End ECM");
}
