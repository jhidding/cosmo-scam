(import (rnrs (6))
        (rnrs io ports (6))
	(scam)
	(cairo)
	(scam ply))

(let* ((T  (vector->list (read-ply-polygon-mesh (cadr (command-line)))))
       (M1 (make-material-linefill-fn
	     (lambda (z normal set-colour set-lw fill stroke)
	       (let* ((c (colour-hsv-gradient (make-colour 'hsva 0.6 1.0 0.2 1.0)
				     	      (make-colour 'hsva 0.2 0.5 1.0 1.0)))
		      (r (abs (a-dot normal (:> 0 0 1)))))
		 (set-colour (c r)) (fill) (set-lw 0.0001)
		 (stroke)))))

       (P (map ($ polygon-add-material -- M1) T))

       (C (camera-transform (:.  -1.0 0.5 0.5)  	 ; position
			    (:.  0 0 0)	 ; target
			    (:>  0 0 1)  	         ; shub
			    parallel-projection))

       (L 600)
       (R (make-svg-renderer C L L "blaauw.svg")))

  (render-do R (lambda (cr)
    (print (cairo-get-line-cap cr) "\n")
    (cairo-set-line-join cr 'round)))
  (render-scale R (* 1/10 L) (* -1/10 L))
  (render-translate R 5 -7)
  (render-scene R P)
  (render-save-png R "blaauw.png")
  (render-finish R))

