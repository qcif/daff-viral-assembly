(function() {
    const KINGDOM_COLORS = {
        'plant': '#4CAF50',
        'animal': '#30b3c7ff',
        'bacteria': '#aa3bb2ff',
        'fungi': '#6f491eff',
        'virus': '#F44336',
        'other': '#9E9E9E'
    };

    function createOverviewPieChart() {
        const kingdoms = [
            { key: 'plant', name: 'Plant' },
            { key: 'animal', name: 'Animal' },
            { key: 'bacteria', name: 'Bacteria' },
            { key: 'fungi', name: 'Fungi' },
            { key: 'virus', name: 'Virus' },
            { key: 'other', name: 'Other' }
        ];

        const labels = [];
        const values = [];
        const colors = [];
        const kingdomKeys = [];

        kingdoms.forEach(kingdom => {
            const kingdomData = krakenTaxaByKingdom[kingdom.key];
            if (kingdomData && kingdomData.read_count > 0) {
                labels.push(kingdom.name);
                values.push(kingdomData.read_count);
                colors.push(KINGDOM_COLORS[kingdom.key]);
                kingdomKeys.push(kingdom.key);
            }
        });

        if (values.length === 0) {
            document.getElementById('kraken-taxa-overview-pie').innerHTML =
                '<div class="alert alert-info">No taxa data available</div>';
            return;
        }

        const data = [{
            labels: labels,
            values: values,
            type: 'pie',
            textinfo: 'label',
            textposition: 'auto',
            hovertemplate: '<b>%{label}</b><br>' +
                          'Read Count: %{value:,}<br>' +
                          'Percentage: %{percent}<br>' +
                          '<extra></extra>',
            marker: {
                colors: colors,
                line: {
                    color: 'white',
                    width: 2
                }
            },
            rotation: -20,
        }];

        const layout = {
            title: {
                text: 'Taxa',
                font: {
                    size: 18,
                    weight: 'bold'
                }
            },
            height: 500,
            width: 500,
            showlegend: false,
            automargin: true
        };

        const config = {
            responsive: true,
            displayModeBar: true,
            displaylogo: false,
            modeBarButtonsToRemove: ['lasso2d', 'select2d', 'toImage']
        };

        Plotly.newPlot('kraken-taxa-overview-pie', data, layout, config);

        document.getElementById('kraken-taxa-overview-pie').on('plotly_click', function(data) {
            console.log("Clicked overview plot");
            console.log(data);
            const pointIndex = data.points[0].i;
            console.log("Clicked point index:", pointIndex);
            const selectedKingdom = kingdomKeys[pointIndex];
            console.log("Selected kingdom:", selectedKingdom);
            showKingdomChart(selectedKingdom);
        });
    }

    function showKingdomChart(kingdomKey) {
        $('.pie-container').hide();
        $(`#kraken-${kingdomKey}-container`).show();
    }

    function createPieChart(containerId, kingdomData, kingdomName) {
        if (!kingdomData || !kingdomData.species || Object.keys(kingdomData.species).length === 0) {
            const capitalizedName = kingdomName.charAt(0).toUpperCase() + kingdomName.slice(1);
            document.getElementById(containerId).innerHTML =
                '<div class="alert alert-info">No data available for ' + capitalizedName + '</div>';
            return;
        }

        const speciesArray = Object.entries(kingdomData.species).map(([species, data]) => ({
            species: species,
            read_count: data.read_count,
            taxon_count: data.taxon_count
        }));

        speciesArray.sort((a, b) => b.read_count - a.read_count);

        const topSpecies = speciesArray.slice(0, N_SPECIES_TO_DISPLAY);
        const otherSpecies = speciesArray.slice(N_SPECIES_TO_DISPLAY);

        const labels = topSpecies.map(item => item.species);
        const values = topSpecies.map(item => item.read_count);

        if (otherSpecies.length > 0) {
            const otherReadCount = otherSpecies.reduce((sum, item) => sum + item.read_count, 0);
            labels.push('Other (' + otherSpecies.length + ' species)');
            values.push(otherReadCount);
        }

        const data = [{
            labels: labels,
            values: values,
            type: 'pie',
            textinfo: null,
            textposition: 'auto',
            hovertemplate: '<b>%{label}</b><br>' +
                          'Read Count: %{value:,}<br>' +
                          'Percentage: %{percent}<br>' +
                          '<extra></extra>',
            marker: {
                line: {
                    color: 'white',
                    width: 2
                }
            },
            rotation: -30,
        }];

        const layout = {
            title: {
                text: kingdomName.charAt(0).toUpperCase() + kingdomName.slice(1) + ' species',
                font: {
                    size: 16,
                    weight: 'bold'
                }
            },
            annotations: [{
              text: `${kingdomData.read_count.toLocaleString()} reads<br>${kingdomData.taxon_count.toLocaleString()} taxa`,
              x: 0.5,
              y: 1.15,
              xref: 'paper',
              yref: 'paper',
              showarrow: false,
              font: { size: 14 }
            }],
            height: 500,
            width: 500,
            showlegend: false,
            automargin: true
        };

        const config = {
            responsive: true,
            displayModeBar: true,
            displaylogo: false,
            modeBarButtonsToRemove: ['lasso2d', 'select2d', 'toImage']
        };

        Plotly.newPlot(containerId, data, layout, config);
    }

    // Create the overview pie chart
    createOverviewPieChart();

    // Create individual kingdom pie charts
    const kingdoms = [
        { id: 'kraken-plant-pie', key: 'plant', name: 'Plant' },
        { id: 'kraken-animal-pie', key: 'animal', name: 'Animal' },
        { id: 'kraken-bacteria-pie', key: 'bacteria', name: 'Bacteria' },
        { id: 'kraken-fungi-pie', key: 'fungi', name: 'Fungi' },
        { id: 'kraken-virus-pie', key: 'virus', name: 'Virus' },
        { id: 'kraken-other-pie', key: 'other', name: 'Other' }
    ];

    kingdoms.forEach(kingdom => {
        createPieChart(kingdom.id, krakenTaxaByKingdom[kingdom.key], kingdom.name);
    });
})();
