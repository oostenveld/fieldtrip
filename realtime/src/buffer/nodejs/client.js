var net = require('net');

var client = new net.Socket();

client.connect(1337, '127.0.0.1', function() {
	client.write('Hello, server! Love, Client.');
});

client.on('data', function(data) {
	console.log('data');
	console.log(data);
	client.destroy(); // kill client after server's response
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

client.on('error', function() {
	console.log('error');
});

client.on('lookup', function() {
	console.log('lookup');
});

client.on('timeout', function() {
	console.log('timeout');
});

