/*
 * Objects and functions to perform the logic of the task
 */

// Static names of HTML elements for use in the functions below
var IMGNAME = "img-canvas";
var MOVNAME = "mov-canvas";
var RBUTTON = "rbutton";
var LBUTTON = "lbutton";
var BBUTTON = "bbutton";
var PROGNAME = "progress";
var DEBUGBTN = "debugbutton";


/*
 * The wrapper class that hold & loads trials and saves information
 */
BalTaskLogic = function(condition,ptobj,endfnc,debugging) {
    if(typeof(debugging)==="undefined"){debugging=false;}
    var trfile = "../static/json/Cond"+(Number(condition)+1)+".json";
    var intrfile = "../static/json/IntroStims.json";

    this.pto = ptobj;
    this.end = endfnc;

    var trlistjson;
    $.ajax({
        dataType: "json",
        url: trfile,
        async: false,
        success: function(data) {
            trlistjson = data;
        },
        error: function(req) {
            throw new Error('Failure to load trial list file');
        }
    });
    var introtrlistjson;
    $.ajax({
        dataType: "json",
        url: intrfile,
        async: false,
        success: function(data) {
            introtrlistjson = data;
        },
        error: function(req) {
            throw new Error('Failure to load intro trial list file');
        }
    });

    var badintro = true;
    while (badintro) {
        this.intrlist = shuffle(introtrlistjson);
        if (this.intrlist[0] !== 'Intro2' & this.intrlist[0] !== 'Intro3') badintro = false;
    }

    this.trlist = shuffle(trlistjson);
    this.tridx = 0;

    this.debugging = debugging;
    this.setup();
};

// Collects all of the document elements and sets the stage
BalTaskLogic.prototype.setup = function() {
    this.img = document.getElementById(IMGNAME);
    this.mov = document.getElementById(MOVNAME);
    this.prog = document.getElementById(PROGNAME);
    this.mov.hidden = true;

    //this.rbutton = document.getElementById(RBUTTON);
    //this.bbutton = document.getElementById(BBUTTON);
    //this.lbutton = document.getElementById(LBUTTON);
    this.canclick = true;

    this.curtr = "";
    this.oninstruct = false;
    this.sttime = 0;

    this.debugbutton = document.getElementById(DEBUGBTN);
    if(!this.debugging) {this.debugbutton.style.visibility = 'hidden';}

    var me = this;

    this.clickable = true;
    keylogfn = function(event) {
        if(me.clickable) {
            console.log(event.keycode);
            if (event.key === "m") { // m: right
                me.clickRight(me);
                me.clickable = false;
            } else if (event.key === "c") { // c: left
                me.clickLeft(me);
                me.clickable = false;
            } else if (event.key === "b") { // b: balance
                me.clickBal(me);
                me.clickable = false;
            }
        }
    }
    document.addEventListener('keydown', keylogfn);
    keyupfn = function(event) {
        me.clickable = true;
    }
    document.addEventListener('keyup', keyupfn);
    //this.rbutton.onclick = function(){me.clickRight(me);};
    //this.lbutton.onclick = function(){me.clickLeft(me);};
    //this.bbutton.onclick = function(){me.clickBal(me);};

    var onload = function() {
        me.allowClick(me);
        me.sttime = new Date().getTime();
    };
    this.img.addEventListener('load',onload,false);
    this.mov.addEventListener('loadeddata',onload,false);
};

BalTaskLogic.prototype.allowClick = function(me) {
    me.canclick = true;
};

// Runs the instruction trials
BalTaskLogic.prototype.startInstruction = function() {
    this.tridx = 0;
    this.img.hidden = true;
    this.mov.hidden = false;
    this.oninstruct = true;
    me = this;
    me.loadInstructTrial(me);
};

BalTaskLogic.prototype.endInstruction = function(me) {
    // Intermezzo instructions
    me.pto.showPage("intermezzo.html");
    var btn = document.getElementById("next");
    btn.addEventListener('click',function(){me.startTrueTrials(me);},false);
};

BalTaskLogic.prototype.startTrueTrials = function(me) {
    me.pto.showPage('stage.html');
    me.tridx = 0;
    me.setup();
    me.loadNewTrial();
};

BalTaskLogic.prototype.loadInstructTrial = function(me) {
    me.curtr = me.intrlist[me.tridx];
    var movnm = "/static/Movies/"+me.curtr+".mp4";
    me.mov.src = movnm;
    me.prog.textContent = ("Progress: "+(me.tridx+1)+" / "+(me.intrlist.length));

};

// Resets the trial with
BalTaskLogic.prototype.loadNewTrial = function() {

    // Load new image onto the screen
    this.curtr = this.trlist[this.tridx];
    var imgnm = "/static/images/"+this.curtr+"_01.jpeg";
    this.img.src = imgnm;

    // Update the progress counter
    this.prog.textContent = ("Progress: "+(this.tridx+1)+" / "+(this.trlist.length));

    //this.canclick = true;
};

// Controls what happens on click
BalTaskLogic.prototype.clickLeft = function(me) {
    if (me.canclick) {
        me.registerResponse('L');
        me.canclick= false;
    }
};
BalTaskLogic.prototype.clickRight = function(me) {
    if (me.canclick) {
        me.registerResponse('R');
        me.canclick = false;
    }
};
BalTaskLogic.prototype.clickBal = function(me) {
    if (me.canclick) {
        me.registerResponse('B');
        me.canclick = false;
    }
};

// Takes in a response, saves to psiturk and gets the next trials
BalTaskLogic.prototype.registerResponse = function(resp) {
    var trtime = new Date().getTime() - this.sttime;
    if(this.oninstruct) {
        this.pto.recordTrialData({
            "phase": "INTROTRIALS",
            "trial": this.curtr,
            "trial_idx": this.tridx,
            "response": resp,
            "rt": trtime,
            "is_instruct": this.oninstruct
        });
        this.tridx = this.tridx + 1;
        var nxtfnc;
        me = this;
        if (this.tridx >= this.intrlist.length) {nxtfnc = function() {me.endInstruction(me);};}
        else {nxtfnc = function() {me.loadInstructTrial(me);};}
        this.mov.addEventListener('ended',nxtfnc,false);
        this.mov.play();
    }
    else {
        this.pto.recordTrialData({
            "phase": "TRIALS",
            "trial": this.curtr,
            "trial_idx": this.tridx,
            "response": resp,
            "rt": trtime,
            "is_instruct": this.oninstruct
        });
        this.tridx = this.tridx + 1;
        if (this.tridx >= this.trlist.length) {this.endExp();}
        else {this.loadNewTrial();}
    }
};

BalTaskLogic.prototype.endExp = function() {
    this.pto.saveData();
    this.end();
};
