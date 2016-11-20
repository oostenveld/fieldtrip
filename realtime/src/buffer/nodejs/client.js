var net = require('net');

var client = new net.Socket();

var get_hdr = Buffer.from([0x00, 0x01, 0x02, 0x01, 0x00, 0x00, 0x00, 0x00]);
var get_dat = Buffer.from([0x00, 0x01, 0x02, 0x02, 0x00, 0x00, 0x00, 0x08]); // followed by begsample and endsample
var get_evt = Buffer.from([0x00, 0x01, 0x02, 0x03, 0x00, 0x00, 0x00, 0x00]);

var get_ok  = Buffer.from([0x00, 0x01, 0x02, 0x04]);
var get_err = Buffer.from([0x00, 0x01, 0x02, 0x05]);


ft_hdr = {
    nchans: 0,
    nsamples: 0,
    nevents: 0,
    fsample: 0,
    data_type: 0
} 

var cmdlist = [];

client.connect(1972, '127.0.0.1', function() {
	cmdlist.push('get_hdr');
	client.write(get_hdr);

	var begsample = new Buffer(4);
	begsample.writeUInt32BE(0);
	var endsample = new Buffer(4);
	endsample.writeUInt32BE(2560);

	cmdlist.push('get_dat');
	client.write(Buffer.concat([get_dat, begsample, endsample]));
});

client.on('data', function(data) {
	console.log('data');
	if (data.slice(0,4).equals(get_ok))
		console.log('get_ok');
	else if (data.slice(0,4).equals(get_err))
		console.log('get_err');
	console.log(data);
});

client.on('close', function() {
	console.log('Connection closed');
});

client.on('connect', function() {
	console.log('connect');
});

client.on('drain', function() {
	console.log('drain');
});

client.on('end', function() {
	console.log('end');
});

client.on('error', function(data) {
	console.log('error');
	console.log(data);
});

client.on('lookup', function() {
	console.log('lookup');
});

client.on('timeout', function() {
	console.log('timeout');
});

