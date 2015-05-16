
package MySubs_IPCELL;
reqiure 5.6.0;

=head1 MyUtils

Plotting functions 

=cut

BEGIN {
}

sub RewriteIpcell
{
  my ($Cactn, $ca_index) = @_;
  my $outputFile = "holzapfel_struct_modify_".$
