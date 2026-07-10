# LaTeX2HTML 2021 (Released January 1, 2021)
# Associate labels original text with physical files.


$key = q/files/;
$external_labels{$key} = "$URL/" . q|node17.html|; 
$noresave{$key} = "$nosave";

$key = q/hort/;
$external_labels{$key} = "$URL/" . q|node19.html|; 
$noresave{$key} = "$nosave";

$key = q/spintransfers/;
$external_labels{$key} = "$URL/" . q|node18.html|; 
$noresave{$key} = "$nosave";

1;


# LaTeX2HTML 2021 (Released January 1, 2021)
# labels from external_latex_labels array.


$key = q/files/;
$external_latex_labels{$key} = q|7|; 
$noresave{$key} = "$nosave";

$key = q/hort/;
$external_latex_labels{$key} = q|9|; 
$noresave{$key} = "$nosave";

$key = q/spintransfers/;
$external_latex_labels{$key} = q|8|; 
$noresave{$key} = "$nosave";

1;

