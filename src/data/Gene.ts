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
    const gene = new Gene(genome.substring(this.start - 1, this.end));

    if (this.strand === 1) {
      return gene;
    } else if (this.strand === -1) {
      return gene.reverseComplement();
    } else {
      throw new Error("Invalid strand");
    }
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

    /*
    Alternative Initiation Codons:

    Trypanosoma: UUA, UUG, CUG
    Leishmania: AUU, AUA
    Tertrahymena: AUU, AUA, AUG
    Paramecium: AUU, AUA, AUG, AUC, GUG, GUA(?)
    */

    const codonTable = {
      4: {
        "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
        "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
        "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "UGU": "C", "UGC": "C", "UGA": "W", "UGG": "W",
        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
      }
    };


    const codons = this.sequence.match(/.{1,3}/g);

    // Check for alternative initiation codons
    const initiationCodons = new Set(["UUA", "UUG", "CUG", "AUU", "AUA", "AUG", "AUC", "GUG", "GUA"]);
    if (initiationCodons.has(codons[0])) {
      codons[0] = "AUG";
    }

    const protein = codons.map(codon => codonTable[table][codon]).join("");

    return new Gene(protein);
  }

  reverseComplement(): Gene {
    const complement = this.sequence.replace(/A/g, "t").replace(/T/g, "a").replace(/C/g, "g").replace(/G/g, "c");
    return new Gene(complement.toUpperCase().split("").reverse().join(""));
  }
}