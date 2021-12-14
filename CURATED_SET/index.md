<script src="msa.min.gz.js"></script>
# Interactive preview of CURATED SET of histone sequences used for HistoneDB 3.0



<div id="yourDiv" class="biojs_msa_div"></div>



<script>
//var msa = require("msa");
        // this is a way how you use a bundled file parser
        // set your custom properties
        // @see: https://github.com/greenify/biojs-vis-msa/tree/master/src/g 
        var opts = {
			el: document.getElementById("yourDiv"),
        	vis: {
            	conserv: false,
            	overviewbox: false,
            	seqlogo: true
        	},
        	conf: {
            	dropImport: true
        	},
        	zoomer: {
        	    menuFontsize: "12px",
        	    autoResize: true
        	}
		};

        // init msa
        var m = new msa.msa(opts);

        var url = "../static/browse/seeds/H3/H3.3.fasta";

        m.u.file.importURL(url, renderMSA);

        function renderMSA() {

            // the menu is independent to the MSA container
            var menuOpts = {
				el: document.getElementById('div'),
				msa: m
			};
            var defMenu = new msa.menu.defaultmenu(menuOpts);
            m.addView("menu", defMenu);

            // call render at the end to display the whole MSA
            m.render();
        }
    </script>
