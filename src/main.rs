mod helpers;
mod primes;
use crate::helpers::*;
use gmp::mpz::*;
use std::fs::File;
use std::io::Write;

// nombres à factoriser pour le challenge
const NUMBERS: [&str; 6] = ["1841831424953158080168", "2886596224401889440451451937054631171754929436102021316072", "2169097458428579869823886745671940306859318149890678788224223007012403826873592", "8967790269526338140041219060275147979025793367237814592721572378013944263405018541606516", "26661657967981291178121484753806009059802561273146506955993517919025168634007436791415970468260073", "9897165300565109444717013704768276163428305757217986431986182768811508720988735974257613403682584583472671"];
const BOUNDS1: [u64; 6] = [
    100000,
    100000,
    3 * 1000000,
    3 * 1000000,
    11 * 1000000,
    43 * 1000000,
];
const BOUNDS2: [u64; 6] = [
    100000,
    100000,
    4592487916,
    4592487916,
    30114149530,
    198654756318,
];

fn main() {
    let mut file = File::create("output.txt").expect("Erreur création fichier");
    for i in 0..5 {
        let n = Mpz::from_str_radix(NUMBERS[i], 10).expect("Error conversion\n");
        let bound1 = Mpz::from(BOUNDS1[i]);
        let bound2 = Mpz::from(BOUNDS2[i]);
        let buf = format!("Integer to factorize : {}\n", n);
        file.write(buf.as_bytes()).expect("Erreur write");
        // println!("Integer to factorize : {}", n);
        factorization(&n, &bound1, &bound2, &mut file);
        file.write(b"\n\nNext Integer\n").expect("Erreur write");
        // print!("\n\nNext Integer\n");
    }
}
