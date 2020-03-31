/*
 * Requires:
 *     psiturk.js
 *     utils.js
 */

var DEBUGGING = false;

// Initalize psiturk object
var psiTurk = new PsiTurk(uniqueId, adServerLoc, mode);

var mycondition = condition;  // these two variables are passed by the psiturk server process
var mycounterbalance = counterbalance;  // they tell you which condition you have been assigned to
// they are not used in the stroop code but may be useful to you

// All pages to be loaded
var pages = [
	"instructions/instruct-1.html",
        "instructions/instruct-wb1.html",
        "instructions/instruct-wb2.html",
        "instructions/instruct-bi1.html",
        "instructions/instruct-bi2.html",
        "instructions/instruct-wi1.html",
        "instructions/instruct-wi2.html",
        "instructions/instruct-shapes.html",
        "instructions/instruct-ssize.html",
        "instructions/instruct-spos.html",
	"instructions/instruct-ready.html",
        "intermezzo.html",
	"stage.html",
	"postquestionnaire.html"
];

psiTurk.preloadPages(pages);

var instructionPages = [ // add as a list as many pages as you like
	"instructions/instruct-1.html",
        "instructions/instruct-wb1.html",
        "instructions/instruct-wb2.html",
        "instructions/instruct-bi1.html",
        "instructions/instruct-bi2.html",
        "instructions/instruct-wi1.html",
        "instructions/instruct-wi2.html",
        "instructions/instruct-shapes.html",
        "instructions/instruct-ssize.html",
        "instructions/instruct-spos.html",
	"instructions/instruct-ready.html"
];

//var instructionPages = ["instructions/instruct-ready.html"]; // For getting past the instructions quickly


/********************
* HTML manipulation
*
* All HTML files in the templates directory are requested
* from the server when the PsiTurk object is created above. We
* need code to get those pages from the PsiTurk object and
* insert them into the document.
*
********************/


/****************
* Questionnaire *
****************/

var Questionnaire = function() {

	var error_message = "<h1>Oops!</h1><p>Something went wrong submitting your HIT. This might happen if you lose your internet connection. Press the button to resubmit.</p><button id='resubmit'>Resubmit</button>";

	record_responses = function() {


		$('textarea').each( function(i, val) {
			psiTurk.recordUnstructuredData(this.id, this.value);
		});
		$('select').each( function(i, val) {
			psiTurk.recordUnstructuredData(this.id, this.value);
		});

	};

	prompt_resubmit = function() {
		replaceBody(error_message);
		$("#resubmit").click(resubmit);
	};

	resubmit = function() {
		replaceBody("<h1>Trying to resubmit...</h1>");
		reprompt = setTimeout(prompt_resubmit, 10000);

		psiTurk.saveData({
			success: function() {
			    clearInterval(reprompt);
                            finish();
			},
			error: prompt_resubmit
		});
	};

	// Load the questionnaire snippet
	psiTurk.showPage('postquestionnaire.html');

	$("#next").click(function () {
	    record_responses();
	    psiTurk.saveData({
            success: function(){
                psiTurk.completeHIT();
            },
            error: prompt_resubmit});
	});


};

var startEnd = function() {
    var Q = new Questionnaire();
};

// Task object to keep track of the current phase
var currentview;

/*******************
 * Run Task
 ******************/
$(window).load( function(){

    psiTurk.recordUnstructuredData("Condition",mycondition);

    var btl;
    var loadExp = function() {
        psiTurk.showPage("stage.html");
        btl = new BalTaskLogic(mycondition,psiTurk,startEnd);
        var dbgbtn = document.getElementById('debugbutton');
        if(!DEBUGGING) {
            dbgbtn.style.visibility = 'hidden';
        } else{
            function dbg() {
                var tnm = btl.curtr;
                var falls = 'fall unknown';
                var beg = tnm.substr(0,2);
                if (beg === 'Ba' | beg === 'CB') {
                    falls = 'balance';
                } else {
                    if (tnm.substr(-1) === 'R') {
                        falls = 'fall right';
                    } else {
                        falls = 'fall left';
                    }
                }
                var alt = "The trial name is " + tnm;
                var alt = alt + "\n \n The beam will " + falls;
                alert(alt);
            };
            dbgbtn.onclick = dbg;
        }
        console.log(btl);
        //btl.startTrueTrials(btl);
        btl.startInstruction();

    };
    psiTurk.doInstructions(instructionPages,loadExp);
	//loadExp();
    /*
    psiTurk.doInstructions(
    	instructionPages, // a list of pages you want to display in sequence
    	function() { currentview = new StroopExperiment(); } // what you want to do when you are done with instructions
    );
    */
});
