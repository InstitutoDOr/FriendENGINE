    #toplevel . -visual truecolor
    wm title . "Model"
    wm iconname . "Model"

    frame .f
    pack .f -expand yes -fill both -in . -side top
    canvas .f.viewport -yscrollcommand ".f.ysbar set" -xscrollcommand ".xsbar set"
    scrollbar .xsbar -command ".f.viewport xview" -orient horizontal
    scrollbar .f.ysbar -command ".f.viewport yview" -orient vertical
    frame .f.viewport.f
    .f.viewport create window 0 0 -anchor nw -window .f.viewport.f
    pack .f.viewport -side left -fill both -expand true -in .f

    set img [ image create photo -file $argv ]
    button .f.viewport.f.btn -command "destroy ." -image $img -borderwidth 0
    pack .f.viewport.f.btn -in .f.viewport.f
    wm geometry . [ image width $img ]x[ image height $img ]