<script src="../static/browse/js/msa%400.4.6.js"></script>
# Interactive preview of CURATED SET of histone sequences used for HistoneDB 3.0


<script>
var msa = require("msa");
// your fasta file (advice: save it in a DOM node)
var fasta = ">seq1\n\
ACTG\n\
>seq2\n\
ACGG\n";

// parsed array of the sequences
var seqs =  msa.io.fasta.parse(fasta);

var m = msa({
     el: rootDiv,
     seqs: seqs
});
m.render();
</script>
