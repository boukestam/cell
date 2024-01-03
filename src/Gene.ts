export interface GenebankFeature {
  name: string;
  start: number;
  end: number;
  strand: 1 | -1;
  type: string;
  notes: {
    [key: string]: string[];
  };
}

export class GeneLocation {
  strand: number;
  start: number;
  end: number;

  constructor(strandOrFeature: number | GenebankFeature, start?: number, end?: number) {
    if (typeof strandOrFeature === 'number') {
      this.strand = strandOrFeature;
      this.start = start!;
      this.end = end!;
    } else {
      const feature = strandOrFeature;
      this.strand = feature.strand;
      this.start = feature.start;
      this.end = feature.end;
    }
  }

  extract(genome: string): Gene {
    return new Gene(genome.substring(this.start, this.end));
  }
}

export class Gene {
  sequence: string;
  length: number;

  constructor(sequence: string) {
    this.sequence = sequence.toUpperCase();
    this.length = sequence.length;
  }

  slice(start: number, end?: number): Gene {
    return new Gene(this.sequence.slice(start, end));
  }

  transcribe(): Gene {
    return new Gene(this.sequence.replace(/T/g, "U"));
  }

  translate(table: number): Gene {
    // Using translation table 4 from NCBI: "Bacterial, Archaeal and Plant Plastid Code"
    // https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG4

    /*
    TTT F Phe      TCT S Ser      TAT Y Tyr      TGT C Cys  
    TTC F Phe      TCC S Ser      TAC Y Tyr      TGC C Cys  
    TTA L Leu i    TCA S Ser      TAA * Ter      TGA W Trp  
    TTG L Leu i    TCG S Ser      TAG * Ter      TGG W Trp  

    CTT L Leu      CCT P Pro      CAT H His      CGT R Arg  
    CTC L Leu      CCC P Pro      CAC H His      CGC R Arg  
    CTA L Leu      CCA P Pro      CAA Q Gln      CGA R Arg  
    CTG L Leu i    CCG P Pro      CAG Q Gln      CGG R Arg  

    ATT I Ile i    ACT T Thr      AAT N Asn      AGT S Ser  
    ATC I Ile i    ACC T Thr      AAC N Asn      AGC S Ser  
    ATA I Ile i    ACA T Thr      AAA K Lys      AGA R Arg  
    ATG M Met i    ACG T Thr      AAG K Lys      AGG R Arg  

    GTT V Val      GCT A Ala      GAT D Asp      GGT G Gly  
    GTC V Val      GCC A Ala      GAC D Asp      GGC G Gly  
    GTA V Val      GCA A Ala      GAA E Glu      GGA G Gly  
    GTG V Val i    GCG A Ala      GAG E Glu      GGG G Gly  
    */

    const codonTable = {
      4: {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "TGT": "C", "TGC": "C", "TGA": "W", "TGG": "W",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
      }
    };

    const codons = this.sequence.match(/.{1,3}/g);
    const protein = codons.map(codon => codonTable[table][codon]).join("");

    return new Gene(protein);
  }

  reverseComplement(): Gene {
    const complement = this.sequence.replace(/A/g, "t").replace(/T/g, "a").replace(/C/g, "g").replace(/G/g, "c");
    return new Gene(complement.toUpperCase().split("").reverse().join(""));
  }
}