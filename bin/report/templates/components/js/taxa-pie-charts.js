function initTaxaPieCharts(prefix, taxaByKingdom, tableSelector) {
    const KINGDOM_COLORS = {
        'plant': '#4CAF50',
        'animal': '#30b3c7ff',
        'bacteria': '#aa3bb2ff',
        'fungi': '#6f491eff',
        'virus': '#F44336',
        'other': '#9E9E9E'
    };

    const KINGDOM_FILTERS = {
        // Define table filters for each kingdom
        // Column indices: 5=full_lineage, 6=entity, 7=domain, 8=kingdom
        'plant': [
            {
                columnIx: 8,
                value: 'viridiplantae'
            },
        ],
        'animal': [
            {
                columnIx: 8,
                value: 'metazoa'
            },
        ],
        'bacteria':  [
            {
                columnIx: 7,
                value: 'bacteria'
            },
        ],
        'fungi': [
            {
                columnIx: 8,
                value: 'fungi'
            },
        ],
        'virus': [
            {
                columnIx: 5,
                value: 'viruses'
            },
        ],
        'other': []
    };

    const rankBtnClass = `${prefix}-rank-btn`;

    // Track current state
    let currentKingdom = null;
    let currentRank = 'family';

    function createOverviewPieChart() {
        const kingdoms = [
            { key: 'plant', name: 'Plant' },
            { key: 'animal', name: 'Animal' },
            { key: 'bacteria', name: 'Bacteria' },
            { key: 'fungi', name: 'Fungi' },
            { key: 'virus', name: 'Virus' },
            { key: 'other', name: 'Other' },
        ];

        const labels = [];
        const values = [];
        const colors = [];
        const kingdomKeys = [];

        kingdoms.forEach(kingdom => {
            const kingdomData = taxaByKingdom[kingdom.key];
            if (kingdomData && kingdomData.read_count > 0) {
                labels.push(kingdom.name);
                values.push(kingdomData.read_count);
                colors.push(KINGDOM_COLORS[kingdom.key]);
                kingdomKeys.push(kingdom.key);
            }
        });

        if (values.length === 0) {
            document.getElementById(`${prefix}-taxa-overview-pie`).innerHTML =
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
                    width: 0,
                }
            },
            rotation: -20,
        }];

        const layout = {
            title: {
                text: 'Classified taxa',
                font: {
                    size: 18,
                    weight: 'bold'
                }
            },
            height: 500,
            width: 500,
            showlegend: false,
            automargin: true,
            margin: { t: 100, b: 100 }
        };

        const config = {
            responsive: true,
            displayModeBar: true,
            displaylogo: false,
            modeBarButtonsToRemove: ['lasso2d', 'select2d', 'toImage']
        };

        Plotly.newPlot(`${prefix}-taxa-overview-pie`, data, layout, config);

        document.getElementById(`${prefix}-taxa-overview-pie`).on('plotly_click', function(data) {
            const pointIndex = data.points[0].i;
            const selectedKingdom = kingdomKeys[pointIndex];
            showKingdomChart(selectedKingdom);
            KINGDOM_FILTERS[selectedKingdom].forEach(
                f => filterTable(f.value, f.columnIx)
            );
        });
    }

    function showKingdomChart(kingdomKey, rank = null) {
        currentKingdom = kingdomKey;
        currentRank = rank || currentRank;

        // Hide default message and show kingdom chart wrapper
        $(`#${prefix}-default-message`).hide();
        $(`#${prefix}-kingdom-chart-wrapper`).show();

        // Hide all kingdom containers, then show the selected one
        $(`#${prefix}-kingdom-chart-container .pie-container`).hide();
        $(`#${prefix}-${kingdomKey}-container`).show();

        updatePieChart(kingdomKey, rank);
    }

    function updatePieChart(kingdomKey, rank = null) {
        rank = rank || currentRank;
        const containerId = `${prefix}-${kingdomKey}-pie`;
        const kingdomData = taxaByKingdom[kingdomKey];
        const capitalizedKingdom = kingdomKey.charAt(0).toUpperCase() + kingdomKey.slice(1);

        if (!kingdomData || !kingdomData[rank] || Object.keys(kingdomData[rank]).length === 0) {
            document.getElementById(containerId).innerHTML =
                '<div class="alert alert-info">No ' + rank + ' data available for ' + capitalizedKingdom + '</div>';
            return;
        }

        const taxaArray = Object.entries(kingdomData[rank]).map(([taxon, data]) => ({
            taxon: taxon,
            read_count: data.read_count,
            taxon_count: data.taxon_count
        }));

        taxaArray.sort((a, b) => b.read_count - a.read_count);
        const labels = taxaArray.map(item => item.taxon);
        const values = taxaArray.map(item => item.read_count);

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
                    width: 0,
                }
            },
            rotation: -30,
        }];

        const layout = {
            title: {
                text: capitalizedKingdom + ' ' + rank,
                font: {
                    size: 16,
                    weight: 'bold'
                }
            },
            annotations: [{
              text: `${kingdomData.read_count.toLocaleString()} reads<br>${kingdomData.taxon_count.toLocaleString()} taxa`,
              x: 0.5,
              y: 1.2,
              xref: 'paper',
              yref: 'paper',
              showarrow: false,
              font: { size: 14 }
            }],
            height: 500,
            width: 500,
            showlegend: false,
            automargin: true,
            margin: { t: 115, b: 115 }
        };

        const config = {
            responsive: true,
            displayModeBar: true,
            displaylogo: false,
            modeBarButtonsToRemove: ['lasso2d', 'select2d', 'toImage']
        };

        Plotly.newPlot(containerId, data, layout, config);

        document.getElementById(containerId).on('plotly_click', function(data) {
            const selectedTaxon = data.points[0].label;
            filterTable(selectedTaxon, 5);
        });
    }

    function filterTable(term, columnId) {
        const dt = $(tableSelector).DataTable();
        dt.search('').columns().search('').draw();
        if (term === 'other' && currentKingdom) {
            // Re-apply the kingdom filter for "other" taxa
            KINGDOM_FILTERS[currentKingdom].forEach(
                f => dt.column(f.columnIx).search(f.value)
            );
        } else {
            dt.column(columnId).search(term);
        }
        dt.draw();
    }

    // Create the overview pie chart
    createOverviewPieChart();

    // Set up rank toggle button event listeners
    $(`.${rankBtnClass}`).on('click', function() {
        const rank = $(this).data('rank');

        // Update button active state
        $(`.${rankBtnClass}`).removeClass('active');
        $(this).addClass('active');

        // Update current rank and refresh the chart
        currentRank = rank;
        if (currentKingdom) {
            updatePieChart(currentKingdom, rank);
        }
    });
}
