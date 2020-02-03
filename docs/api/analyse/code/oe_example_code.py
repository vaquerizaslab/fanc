# start snippet oe setup
import fanc
import matplotlib.pyplot as plt

hic_500kb = fanc.load("output/hic/binned/fanc_example_500kb.hic")
# end snippet oe setup

hic_500kb.close()
# start snippet oe append
hic_500kb = fanc.load("output/hic/binned/fanc_example_500kb.hic", mode='a')
# end snippet oe append


# start snippet oe basic
intra_expected, intra_expected_chromosome, inter_expected = hic_500kb.expected_values()
# end snippet oe basic

# start snippet oe ddplot
# obtain bin distances
bin_size = hic_500kb.bin_size
distance = list(range(0, bin_size * len(intra_expected_chromosome['chr19']), bin_size))

# plot expected values
fig, ax = plt.subplots()
plt.plot(distance, intra_expected_chromosome['chr19'])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel("Distance")
ax.set_ylabel("Average contacts")
plt.show()
# end snippet oe ddplot
fig.savefig('../fanc-doc/api/analyse/images/oe_500kb.png')


# start snippet oe nonorm
intra_expected_nonorm, intra_expected_chromosome_nonorm, inter_expected_nonorm = hic_500kb.expected_values(norm=False)

# obtain bin distances
bin_size = hic_500kb.bin_size
distance = list(range(0, bin_size * len(intra_expected_nonorm), bin_size))

# plot expected values
fig, ax = plt.subplots()
plt.plot(distance, intra_expected_nonorm)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel("Distance")
ax.set_ylabel("Average contacts (not normalised)")
plt.show()
# end snippet oe nonorm
fig.savefig('../fanc-doc/api/analyse/images/oe_500kb_nonorm.png')
